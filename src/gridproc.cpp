#include"gridproc.h"
#include"omp.h"
#include"mkl.h"
void wrapx(float* x1, float* x2, int J, int Nthreads)
{
	omp_set_num_threads(Nthreads);
#pragma omp parallel for
	for(int i=0;i<J;i++)
	{
		x1[i]=fmod(x1[i]-(int)x1[i]+0.5f+1,1)-0.5f;
		x2[i]=fmod(x2[i]-(int)x2[i]+0.5f+1,1)-0.5f;
	}
}

void quicksort_inds(int* arr, int* inds, int left, int right) 
{
	int i=left,j=right;
	int tmp;
	int pivot = arr[(left+right)/2];
	while (i<=j) {
		while (arr[i] < pivot)i++;
		while (arr[j] > pivot)j--;
		if (i <= j) {
			tmp = arr[i];arr[i] = arr[j];arr[j] = tmp;
			tmp = inds[i];inds[i] = inds[j];inds[j] = tmp;
			i++;j--;
		}
	};
	if (left < j)quicksort_inds(arr,inds,left,j);
	if (i < right)quicksort_inds(arr,inds,i,right);
}

void reorder(float *x,float *y,float *ax,float *ay, float* tmpf,int* inds,int J)
{
	cblas_sgthr(J,y,tmpf,inds);
	cblas_scopy (J,tmpf,1,y,1);
	cblas_sgthr(J,x,tmpf,inds);
	cblas_scopy (J,tmpf,1,x,1);
	cblas_sgthr(J,ax,tmpf,inds);
	cblas_scopy (J,tmpf,1,ax,1);
	cblas_sgthr(J,ay,tmpf,inds);
	cblas_scopy (J,tmpf,1,ay,1);
}
void reorderf(cplxf *f,cplxf* tmpf,int* inds, int J)
{
	cblas_cgthr(J,f,tmpf,inds);
	cblas_ccopy (J,tmpf,1,f,1);
}
void reorderf_back(cplxf *f,cplxf* tmpf,int* inds,int J)
{
	cblas_csctr(J,f,inds,tmpf);
	cblas_ccopy (J,tmpf,1,f,1);
}

void sort_by_key(int2* st_num, float *x,float *y,float *ax,float *ay, int* xi_ind, int J,int Msize)
{	
	for(int i=0;i<Msize;i++){st_num[i].x=0;st_num[i].y=0;}
	int pos=0;int count=0;int st;
	while(pos<J)
	{
		st=xi_ind[pos];
		st_num[st].x=pos;
		while(pos<J&&st==xi_ind[pos])
		{
			count++;
			pos++;
		}
		st_num[st].y=count;
		count=0;
	}
}