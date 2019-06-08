#include"expradon.h"
#include"gridproc.h"
#include<string.h>
#include<mex.h>
expRadon::expRadon(int* pints,float* pfloats)
{
	//mexPrintf("init class\n");
	N=pints[0];Nthreads=pints[1];nu=pints[2];
	A=pfloats[0];epsilon=pfloats[1];

	//precompute parameters
	mu=-log((float)epsilon)/(nu*N*N*(nu-1))+(2*nu-1)/(float)(2*nu*N*N*(nu-1))*log((float)A);
	float Te=1/PIf*sqrtf(-mu*log((float)epsilon)+(mu*N)*(mu*N)/4.0f+log((float)A)*log((float)A)/(4.0f*N*N));
	M=(int)ceil(nu*N*Te);

	//strides for fft
	strides[1]=1;strides[2]=2*M+nu*N;strides[0]=M*strides[1]+M*strides[2];

	//cplx array for fft
	G=new cplxf[(2*M+nu*N)*(2*M+nu*N)];memset(G,0,(2*M+nu*N)*(2*M+nu*N)*sizeof(cplxf));

	//plans for fft
	MKL_LONG Nf[2];Nf[0]=nu*N;Nf[1]=nu*N;
    DftiCreateDescriptor(&fwd_plan, DFTI_SINGLE, DFTI_COMPLEX, 2, Nf);
	DftiSetValue(fwd_plan,DFTI_INPUT_STRIDES,strides);
	DftiSetValue(fwd_plan,DFTI_OUTPUT_STRIDES,strides);
	DftiCommitDescriptor(fwd_plan);

	//temporarily arrays
	phiaxt = new cplxf[(2*M+1)*Nthreads];
	phiayt = new cplxf[(2*M+1)*Nthreads];
	Gphiaxt = new cplxf[(2*M+1)*Nthreads];

	init_grids=0;
}
expRadon::~expRadon()
{
//	mexPrintf("delete class\n");
	if(init_grids)
	{
		delete[] x;
		delete[] y;
		delete[] ax;
		delete[] ay;
		delete[] inds;
		delete[] st_num;
		delete[] tmparray;
		delete[] F;
		delete[] weights;
		DftiFreeDescriptor(&fwd_plan1d);
		init_grids=0;
	}
	delete[] G;
	delete[] phiaxt;
	delete[] phiayt;
	delete[] Gphiaxt;   
	DftiFreeDescriptor(&fwd_plan); 
	mkl_free_buffers();
}

void expRadon::fftshift(cplxf* G)
{
	omp_set_num_threads(Nthreads);
#pragma omp parallel for
	for(int j2=0;j2<nu*N;j2++) 
	{
		int k1=(1-2*((j2+1)&1));
		for(int j1=(k1==-1);j1<nu*N;j1+=2)
		{
			G[strides[0]+j1+(j2)*strides[2]].x=-G[strides[0]+j1+(j2)*strides[2]].x; //strides[1]==1;
			G[strides[0]+j1+(j2)*strides[2]].y=-G[strides[0]+j1+(j2)*strides[2]].y;
		}
	}
}
void expRadon::fftshift1d(cplxf* G)
{
	for(int j1=0;j1<Nsf;j1+=2)
	{
		G[j1].x=-G[j1].x;
		G[j1].y=-G[j1].y;
	}
}

void expRadon::expfft1d(float* g,float* f)
{	
	omp_set_num_threads(Nthreads);
	mkl_set_num_threads(1);
#pragma omp parallel for
	for(int i=0;i<Ntheta;i++)
	{
		for(int j=0;j<Nsf;j++)
		{
			tmparray[i*Nsf+j].x=f[i*Nsf+j];
			tmparray[i*Nsf+j].y=0.0;
		}
		fftshift1d((cplxf*)&tmparray[i*Nsf].x);
		mkl_set_num_threads(1);
		DftiComputeForward( fwd_plan1d, &tmparray[i*Nsf]);
		fftshift1d((cplxf*)&tmparray[i*Nsf].x);
		memcpy(&g[2*i*Nsf],&tmparray[i*Nsf].x,Nsf*sizeof(cplxf));
	}
	mkl_set_num_threads(Nthreads);
}
void expRadon::expifft1d(float* g,float* f)
{
	memcpy(tmparray,f,2*Nsf*Ntheta*sizeof(float));
	omp_set_num_threads(Nthreads);
	mkl_set_num_threads(1);
#pragma omp parallel for
	for(int i=0;i<Ntheta;i++)
	{
		fftshift1d((cplxf*)&tmparray[i*Nsf].x);
		DftiComputeBackward( fwd_plan1d, &tmparray[i*Nsf].x);
		fftshift1d((cplxf*)&tmparray[i*Nsf].x);
		for(int j=0;j<Nsf;j++) g[i*Nsf+j]=tmparray[i*Nsf+j].x;//take real
	}
	mkl_set_num_threads(Nthreads);
}

void expRadon::eq2us(float *ut,float *f)
{	
	memset(G,0,(2*M+nu*N)*(2*M+nu*N)*sizeof(cplxf));
	//mul by weight function
	for(int j2=-N/2;j2<N/2;j2++)for(int j1=-N/2;j1<N/2;j1++)
	{
		int indx=(nu*N/2+j1);
		int indy=(nu*N/2+j2);
		G[strides[0]+indx*strides[1]+strides[2]*indy].x=f[(j2+N/2)*N+j1+N/2]*expf(mu*j1*j1+mu*j2*j2)/(nu*N*nu*N);
		G[strides[0]+indx*strides[1]+strides[2]*indy].y=0;
	}
	mkl_set_num_threads(Nthreads);//
	
	fftshift(G);
    DftiComputeForward(fwd_plan, &G[0]);
	fftshift(G);
	mkl_set_num_threads(1);//
	gather((cplxf*)ut,G);
}

void expRadon::us2eq(float *ut,float *fin,float* w,int* ids)
{
	for(int j=0;j<Ns;j++)
	{
		weights[j].x=w[j];
		weights[j].y=0;
	}
	omp_set_num_threads(Nthreads);
	mkl_set_num_threads(1);
#pragma omp parallel for
	for(int i=0;i<Ntheta;i++)
	{
		//copy by index
		cblas_cgthr(Ns,(cplxf*)&fin[i*2*Nsf],&F[i*Ns].x,ids);
		//mul weights
		vcMul(Ns,(MKL_Complex8*)&F[i*Ns],(MKL_Complex8*)weights,(MKL_Complex8*)&F[i*Ns]);
	}


	memset(G,0,(2*M+nu*N)*(2*M+nu*N)*sizeof(cplxf));
	scatter(G,(cplxf*)F);

	mkl_set_num_threads(Nthreads);
	fftshift(G);
	DftiComputeForward(fwd_plan, &G[0]);
	fftshift(G);
	mkl_set_num_threads(1);//	
#pragma omp parallel for
	for(int j2=-N/2;j2<N/2;j2++)for(int j1=-N/2;j1<N/2;j1++)
	{
		int indx=(nu*N/2+j1);
		int indy=(nu*N/2+j2);
		ut[(j2+N/2)*N+j1+N/2]=G[strides[0]+indx*strides[1]+strides[2]*indy].x*expf(mu*j1*j1+mu*j2*j2)/(nu*N*nu*N);
	}
}



void expRadon::set_grids(float* x_,float* y_,float* ax_,float* ay_,int* pints)
{
	if(init_grids)
	{
		delete[] x;
		delete[] y;
		delete[] ax;
		delete[] ay;
		delete[] inds;
		delete[] st_num;
		delete[] tmparray;
		delete[] F;
		delete[] weights;
		DftiFreeDescriptor(&fwd_plan1d);
		init_grids=0;
	}
	Ns=pints[0];Ntheta=pints[1];Nsf=pints[2];
	//copy grids
	x=new float[Ns*Ntheta];memcpy(x,x_,Ns*Ntheta*sizeof(float));
	y=new float[Ns*Ntheta];memcpy(y,y_,Ns*Ntheta*sizeof(float));
	ax=new float[Ns*Ntheta];memcpy(ax,ax_,Ns*Ntheta*sizeof(float));
	ay=new float[Ns*Ntheta];memcpy(ay,ay_,Ns*Ntheta*sizeof(float));
	inds=new int[Ns*Ntheta];
	st_num= new int2[nu*nu*N*N];
	tmparray = new cplxf[Ntheta*Nsf];
	F=new cplxf[Ns*Ntheta];
	weights=new cplxf[Ns];

	//wrap values out of the interval [-0.5 0.5)
	wrapx(x,y,Ntheta*Ns,Nthreads);

	//sort points
	int* xi_ind=new int[Ns*Ntheta];
	for(int i=0;i<Ns*Ntheta;i++)
	{
		int indx=(int)(nu*N*x[i])+nu*N/2;
		int indy=(int)(nu*N*y[i])+nu*N/2;
		xi_ind[i]=indx+nu*N*indy;//nearest int
		inds[i]=i;
	}
	quicksort_inds(xi_ind,inds,0,Ns*Ntheta-1);
	reorder(x,y,ax,ay,(float*)tmparray,inds,Ns*Ntheta);
	sort_by_key(st_num,x,y,ax,ay,xi_ind,Ns*Ntheta,nu*nu*N*N);
	delete[] xi_ind;
	
	//1d ffts
	DftiCreateDescriptor( &fwd_plan1d, DFTI_SINGLE,DFTI_COMPLEX, 1, Nsf);
	DftiSetValue(fwd_plan1d,DFTI_BACKWARD_SCALE,1.0/(float)Nsf);
	DftiSetValue(fwd_plan1d,DFTI_FORWARD_SCALE,1.0/(float)Nsf);
	DftiCommitDescriptor(fwd_plan1d);
	
	init_grids=1;
}
void expRadon::getSizes(size_t* sizes)
{
	sizes[0]=N;
	sizes[1]=N;
	sizes[2]=Ns;
	sizes[3]=Ntheta;
	sizes[4]=2*Ns;
	sizes[5]=Ntheta;
	sizes[6]=Nsf;
	sizes[7]=Ntheta;
	sizes[8]=2*Nsf;
	sizes[9]=Ntheta;
}