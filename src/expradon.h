#pragma once
#include"mkl.h"
#include"omp.h"
#include<mex.h>
#include"config.h"
class expRadon
{
	int N;int nu;
	int Ntheta;	int Ns;	int Nsf;
	int Nthreads;	
	int M;
	MKL_LONG strides[3];
	float A; float epsilon;
	float mu;
	float* x; float* y;
	float* ax; float* ay;
	cplxf* G;
	cplxf* F;
	//fft handles
	DFTI_DESCRIPTOR_HANDLE fwd_plan;
	DFTI_DESCRIPTOR_HANDLE fwd_plan1d;

	int* inds;cplxf* tmparray;int2* st_num;
	bool init_grids; //if grids are initialized
	
	cplxf* phiaxt; cplxf* phiayt;cplxf* Gphiaxt;
	
	cplxf* weights;

public:
	expRadon(int* pints,float* pfloat);
	~expRadon(void);
	void fftshift(cplxf* G);
	void fftshift1d(cplxf* G);
	void gather(cplxf* ut,cplxf* G);
	void scatter(cplxf* G,cplxf* f);
	void eq2us(float* ut,float* f);//float instead of cplxf for matlab interface compatibility
	void us2eq(float *ut,float *f,float* weights, int* ids);
	void expfft1d(float* g,float* f);
	void expifft1d(float* g,float* f);
	void getSizes(size_t* sizes);
	void printParams();	
	void set_grids(float* x_,float* y_,float* ax_,float* ay_,int* pints);
};

