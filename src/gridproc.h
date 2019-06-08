#include"config.h"
void wrapx(float* x1, float* x2,int J, int Nthreads);
void quicksort_inds(int* arr, int* inds, int left, int right);
void reorder(float *x,float *y,float *ax,float *ay, float *tmpf,int* inds,int J);
void reorderf(cplxf *f,cplxf* tmpf,int* inds,int J);
void reorderf_back(cplxf *f,cplxf* tmpf,int* inds,int J);
void sort_by_key(int2* st_num, float *x,float *y,float *ax,float *ay, int* xi_ind, int J,int Msize);