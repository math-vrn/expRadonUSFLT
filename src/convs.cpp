#include"expradon.h"
#include"gridproc.h"

void expRadon::gather(cplxf* f,cplxf* G)
{
	//Wrap boundary
#pragma omp parallel for collapse(2)
		for (int j2=-M;j2<nu*N+M;j2++)for (int j1=-M;j1<nu*N+M;j1++)if ((j1<0)||(j1>=nu*N)||(j2<0)||(j2>=nu*N))
		{
			int ind1=strides[0]+strides[1]*j1+strides[2]*j2;
			int ind2=strides[0]+strides[1]*((j1+nu*N)%(nu*N))+strides[2]*((j2+nu*N)%(nu*N));
			G[ind1].x=G[ind2].x;
			G[ind1].y=G[ind2].y;
		}		

		const cplxf alpha={1.0, 0.0};
		const cplxf beta={0.0, 0.0};
		const float pidmu1=PIf/mu;
		const float pidmu2=PIf/mu;
		const float pi2dmu1=PIf*PIf/mu;
		const float pi2dmu2=PIf*PIf/mu;
		const float c21p=expf(-2*pi2dmu1/(float)(nu*N*nu*N)); 
		const float c22p=expf(-2*pi2dmu2/(float)(nu*N*nu*N));

		omp_set_num_threads(Nthreads);
		mkl_set_num_threads(1);
#pragma omp parallel for
		for(int k=0;k<Ns*Ntheta;k++)
		{
			int ithread=omp_get_thread_num();

			cplxf* phiax = (cplxf*)&phiaxt[(2*M+1)*ithread];
			cplxf* phiay = (cplxf*)&phiayt[(2*M+1)*ithread];
			cplxf* Gphiax = (cplxf*)&Gphiaxt[(2*M+1)*ithread];

			f[k].x=0.0;f[k].y=0.0;

			int ell1=(int)(nu*N*x[k]);
			int ell2=(int)(nu*N*y[k]);
			float dax= ax[k]/(2*PIf);
			float day= ay[k]/(2*PIf);
			//phiax

			cplxf pow1;
			cplxf phiaxi;
			cplxf b1;
			cplxf tmp;
			float dx=((ell1-M)/(float)(nu*N)-x[k]);
			pow1.x=-pi2dmu1*(dx*dx-dax*dax);
			pow1.y=pi2dmu1*2*dx*dax;


			phiaxi.x=sqrtf(pidmu1*pidmu2)*expf(pow1.x)*cosf(pow1.y);
			phiaxi.y=sqrtf(pidmu1*pidmu2)*expf(pow1.x)*sinf(pow1.y);
			//pow1.x=-2*pi2dmu1/(nu*N)*( (ell1-M+0.5)/(float)(nu*N) - x[k] );
			pow1.x=-2*pi2dmu1/(nu*N)*( (ell1-M)/(float)(nu*N) - x[k] +0.5f/(float)(nu*N ));
			pow1.y=ax[k]*pidmu1/(nu*N);
			b1.x=expf(pow1.x)*cosf(pow1.y);
			b1.y=expf(pow1.x)*sinf(pow1.y);
			for(int j1=0;j1<2*M+1;j1++)
			{
				tmp=phiaxi;	
				phiax[j1]=phiaxi;
				phiaxi.x=tmp.x*b1.x-tmp.y*b1.y;
				phiaxi.y=tmp.y*b1.x+tmp.x*b1.y;
				b1.x*=c21p;
				b1.y*=c21p;
			}
			//phiay

			cplxf pow2;
			cplxf phiayi;
			cplxf b2;
			float dy=((ell2-M)/(float)(nu*N)-y[k]);
			pow2.x=-pi2dmu2*(dy*dy-day*day);
			pow2.y=pi2dmu2*2*dy*day;

			phiayi.x=expf(pow2.x)*cosf(pow2.y);
			phiayi.y=expf(pow2.x)*sinf(pow2.y);
			//b2
			//pow2.x=-2*pi2dmu2/(nu*N)*( (ell2-M+0.5)/(float)(nu*N) - y[k] );
			pow2.x=-2*pi2dmu2/(nu*N)*( (ell2-M)/(float)(nu*N) - y[k] +0.5f/(float)(nu*N ));
			pow2.y=ay[k]*pidmu2/(nu*N);
			b2.x=expf(pow2.x)*cosf(pow2.y);
			b2.y=expf(pow2.x)*sinf(pow2.y);
			for(int j2=0;j2<2*M+1;j2++)
			{
				tmp=phiayi;	
				phiay[j2]=phiayi;
				phiayi.x=tmp.x*b2.x-tmp.y*b2.y;
				phiayi.y=tmp.y*b2.x+tmp.x*b2.y;
				b2.x*=c22p;
				b2.y*=c22p;
			}
			int indx=nu*N/2+(ell1-M);
			int indy=nu*N/2+(ell2-M);

			cblas_cgemv(CblasRowMajor,CblasNoTrans,2*M+1,2*M+1,&alpha,
				(cplxf*)&G[strides[0]+indx*strides[1]+strides[2]*indy], strides[2],phiax,1,&beta,Gphiax,1);//y := alpha*A*x + beta*y,

			cplxf fp={0,0};
			cblas_cdotu_sub(2*M+1,phiay,1,Gphiax,1,&fp);

			f[k].x+=fp.x;
			f[k].y+=fp.y;

		}
		mkl_set_num_threads(Nthreads);//return parallel MKL
		reorderf_back(f,tmparray,inds,Ns*Ntheta);		
}

void expRadon::scatter(cplxf *G, cplxf *f)
{
	float norm=0;
	mkl_set_num_threads(Nthreads);//return parallel MKL
	reorderf(f,tmparray,inds,Ns*Ntheta);

	int P1gg,P2gg;
	P1gg=64;P2gg=64;
	int Ntp=(int)(ceil(  ceil(nu*N/(float)(2*M+P1gg)) * ceil(nu*N/(float)(2*M+P2gg)) /(float)Nthreads)+0.5);


	int stepx=2*M+P1gg;	
	int stepy=2*M+P2gg;    

	const cplxf alpha={1.0,0.0};
	const float pidmu1=PIf/mu;
	const float pidmu2=PIf/mu;
	const float pi2dmu1=PIf*PIf/mu;
	const float pi2dmu2=PIf*PIf/mu;
	const float c21p=expf(-2*pi2dmu1/(float)(nu*N*nu*N)); 
	const float c22p=expf(-2*pi2dmu2/(float)(nu*N*nu*N));	
	double t_start = omp_get_wtime();
	omp_set_num_threads(Nthreads);
	mkl_set_num_threads(1);
#pragma omp parallel
	{
		for(int istepy=0;istepy<stepy;istepy+=P2gg)
			for(int istepx=0;istepx<stepx;istepx+=P1gg)
			{
				int P1g=P1gg;
				int P2g=P2gg;
				if(istepx+P1gg>=stepx) P1g=stepx-istepx;
				if(istepy+P2gg>=stepy) P2g=stepy-istepy;

				int ithread=omp_get_thread_num();
				cplxf* phiax = (cplxf*)&phiaxt[(2*M+1)*ithread];
				cplxf* phiay = (cplxf*)&phiayt[(2*M+1)*ithread];
				cplxf* Gphiax = (cplxf*)&Gphiaxt[(2*M+1)*ithread];
				cplxf pow1;
				cplxf phiaxi;
				cplxf b1;
				cplxf tmp;
				cplxf pow2;
				cplxf phiayi;
				cplxf b2;

#pragma omp for collapse(2) schedule(dynamic)
				for(int j=istepy;j<nu*N;j+=stepy)for(int i=istepx;i<nu*N;i+=stepx)
				{
					int P1=P1g;
					int P2=P2g;
					if(i+P1g>=nu*N) P1=nu*N-i;
					if(j+P2g>=nu*N) P2=nu*N-j;

					int indgg=strides[0]+(j-M)*strides[2]+i-M;

					for(int p2=0;p2<P2;p2++)
					{
						int indg=indgg+p2*strides[2];
						for(int p1=0;p1<P1;p1++)
						{

							int ell1=i-nu*N/2+p1;
							int ell2=j-nu*N/2+p2;//fixed
							cplxf* Gl=(cplxf*)&G[indg+p1];
							int2 st_numl=st_num[(j+p2)*(nu*N)+i+p1];
							
							for(int indk=0;indk<st_numl.y;indk++)
							{
								int k=st_numl.x+indk;
								float dax= ax[k]/(2*PIf);
								float day= ay[k]/(2*PIf);
								//phiax
								float dx=((ell1-M)/(float)(nu*N)-x[k]);
								pow1.x=-pi2dmu1*(dx*dx-dax*dax);
								pow1.y=pi2dmu1*2*dx*dax;
								phiaxi.x=sqrtf(pidmu1*pidmu2)*expf(pow1.x)*cosf(pow1.y);
								phiaxi.y=sqrtf(pidmu1*pidmu2)*expf(pow1.x)*sinf(pow1.y);
								pow1.x=-2*pi2dmu1/(nu*N)*( (ell1-M)/(float)(nu*N) - x[k] +0.5f/(float)(nu*N ));
								pow1.y=ax[k]*pidmu1/(nu*N);
								b1.x=expf(pow1.x)*cosf(pow1.y);
								b1.y=expf(pow1.x)*sinf(pow1.y);
								for(int j1=0;j1<2*M+1;j1++)
								{
									tmp=phiaxi;	
									phiax[j1]=phiaxi;
									phiaxi.x=tmp.x*b1.x-tmp.y*b1.y;
									phiaxi.y=tmp.y*b1.x+tmp.x*b1.y;
									b1.x*=c21p;
									b1.y*=c21p;
								}
								//phiay

								float dy=((ell2-M)/(float)(nu*N)-y[k]);
								pow2.x=-pi2dmu2*(dy*dy-day*day);
								pow2.y=pi2dmu2*2*dy*day;
								phiayi.x=expf(pow2.x)*cosf(pow2.y);
								phiayi.y=expf(pow2.x)*sinf(pow2.y);
								//b2
								pow2.x=-2*pi2dmu2/(nu*N)*( (ell2-M)/(float)(nu*N) - y[k] +0.5f/(float)(nu*N ));
								pow2.y=ay[k]*pidmu2/(nu*N);
								b2.x=expf(pow2.x)*cosf(pow2.y);
								b2.y=expf(pow2.x)*sinf(pow2.y);
								for(int j2=0;j2<2*M+1;j2++)
								{
									tmp=phiayi;	
									phiay[j2]=phiayi;
									phiayi.x=tmp.x*b2.x-tmp.y*b2.y;
									phiayi.y=tmp.y*b2.x+tmp.x*b2.y;
									b2.x*=c22p;
									b2.y*=c22p;
								}

								cblas_cscal(2*M+1,&f[k].x,phiax,1);

								cblas_cgeru(CblasRowMajor,2*M+1,2*M+1,&alpha,phiay,1,phiax,1,Gl,strides[2]);
							}
						}
					}
				}
			}
	}

	//Wrap boundary
#pragma omp parallel for collapse(2)
	for (int j2=-M;j2<nu*N+M;j2++)for (int j1=-M;j1<nu*N+M;j1++)if ((j1<0)||(j1>=nu*N)||(j2<0)||(j2>=nu*N))
	{
		int ind1=strides[0]+strides[1]*j1+strides[2]*j2;
		int ind2=strides[0]+strides[1]*((j1+nu*N)%(nu*N))+strides[2]*((j2+nu*N)%(nu*N));
		G[ind2].x+=G[ind1].x;
		G[ind2].y+=G[ind1].y;
	}
	mkl_set_num_threads(Nthreads);
}
