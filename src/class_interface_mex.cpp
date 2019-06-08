#include "mex.h"
#include "class_handle.hpp"
#include "expradon.h"
// The class that we are interfacing to
#include<iostream>
using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    if (nrhs == 0) {
        mexUnlock();
        return;
    }
	// Get the command string
	char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

	// New
	if (!strcmp("new", cmd)) {
		// Check parameters
		if (nlhs != 1 || nrhs!=3)
			mexErrMsgTxt("New: Unexpected arguments.");
		// Return a handle to a new C++ instance
        int *a = (int*)mxGetData(prhs[1]);
        float *b = (float*)mxGetData(prhs[2]);        
		expRadon* expRadon0=new expRadon((int*)mxGetData(prhs[1]),(float*)mxGetData(prhs[2]));
		plhs[0] = convertPtr2Mat<expRadon>(expRadon0);

		return;
	}

	// Check there is a second input, which should be the class instance handle
	if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");

	if (!strcmp("delete", cmd)) {
		// Destroy the C++ object
		destroyObject<expRadon>(prhs[1]);
		if (nlhs != 0 || nrhs != 2)
			mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
		return;
	}

	// Get the class instance pointer from the second input
	expRadon *expRadon0 = convertMat2Ptr<expRadon>(prhs[1]);
	if (!strcmp("eq2us", cmd)) {//eq2us(cplxf *ut,cplxf *f)
		if (nlhs > 1 || nrhs != 3)
			mexErrMsgTxt("eq2us: Unexpected arguments.");
		size_t sizes[10];;
		expRadon0->getSizes(sizes);
		plhs[0] = mxCreateNumericArray(2, &sizes[4], mxSINGLE_CLASS, mxREAL);
		expRadon0->eq2us((float*)mxGetData(plhs[0]),(float*)mxGetData(prhs[2]));
		return;
	}
	if (!strcmp("us2eq", cmd)) {//
		if (nlhs > 1 || nrhs != 5)
			mexErrMsgTxt("us2eq: Unexpected arguments.");
		size_t sizes[10];;
		expRadon0->getSizes(sizes);
		plhs[0] = mxCreateNumericArray(2, &sizes[0], mxSINGLE_CLASS, mxREAL);
		expRadon0->us2eq((float*)mxGetData(plhs[0]),(float*)mxGetData(prhs[2]),(float*)mxGetData(prhs[3]),(int*)mxGetData(prhs[4]));
		return;
	}
	if (!strcmp("expifft1d", cmd)) {
		if (nlhs > 1 || nrhs != 3)
			mexErrMsgTxt("expifft1d: Unexpected arguments.");
		size_t sizes[10];;
		expRadon0->getSizes(sizes);
		plhs[0] = mxCreateNumericArray(2, &sizes[6], mxSINGLE_CLASS, mxREAL);
		expRadon0->expifft1d((float*)mxGetData(plhs[0]),(float*)mxGetData(prhs[2]));
		return;
	}
	if (!strcmp("expfft1d", cmd)) {
		if (nlhs > 1 || nrhs != 3)
			mexErrMsgTxt("expfft1d: Unexpected arguments.");
		size_t sizes[10];;
		expRadon0->getSizes(sizes);
		plhs[0] = mxCreateNumericArray(2, &sizes[8], mxSINGLE_CLASS, mxREAL);
		expRadon0->expfft1d((float*)mxGetData(plhs[0]),(float*)mxGetData(prhs[2]));
		return;
	}
	if (!strcmp("set_grids", cmd)) {
		if (nlhs > 1 || nrhs != 7)
			mexErrMsgTxt("set_grids: Unexpected arguments.");
		expRadon0->set_grids((float*)mxGetData(prhs[2]),(float*)mxGetData(prhs[3]),(float*)mxGetData(prhs[4]),(float*)mxGetData(prhs[5]),(int*)mxGetData(prhs[6]));
		return;
	}

	// Got here, so command not recognized
	mexErrMsgTxt("Command not recognized.");
}
