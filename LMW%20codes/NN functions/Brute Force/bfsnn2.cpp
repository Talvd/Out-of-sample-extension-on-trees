// bfsnn (Brute force nearest neighbor search)
//
//
//
// Each column correponds to a point. Each row is a coordinate.
//
// Yoel Shkolnisky, May 2007.

//#define DEBUG

#include <mex.h>
#include <math.h>
#include <string.h>
#include "QuickSort.h"


#define IMAGE_NAME "bfsnn"

#define for if(0);else for

/*
* MATLAB utlity routines
*/


bool GetInteger(const mxArray *a, int *value)
{
    if( !mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a))
        return false;
    if(mxGetM(a)*mxGetN(a) > 1)
        return false;
    double *pr = mxGetPr(a);
    
    // check to see that the value is actually an integer
    if( floor(pr[0])!=pr[0])
        return false;
    
    *value = (int)pr[0];
    return true;
}


void usage()
{
	printf("Brute-force search of nearest neighbors.\n");
	printf("Find for each point in Y the NN nearest points among the points in X.\n");	
	printf("\n");
	printf("[idx,dists]=bfsnn(X,Y,NN)\n");
	printf("\n");
	printf("   X      M by N dataset. N data points. M coordinates to each data point.\n");
	printf("   Y      M by K dataset. K data points. M coordinates to each data point.\n");
	printf("          If Y is missing then Y=X.\n");
	printf("   NN     Number of nearest neighbors to find for each data point\n");
	printf("\n");
	printf("Returned values:\n");
	printf("   idx    NN by K array. Column i contains the indices of the NN nearest points to i in X.\n");
	printf("   dists  NN by K array. Column i contains the squared distances to NN nearest neighbors in X.\n");
}



void error_msg(const char* msg)
{
	printf("%s: %s",IMAGE_NAME,msg);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	
    /* if there are no input arguments, print out help information */
    if(nrhs==0) {
        usage();
        return;
    }       
    
    if(nrhs < 2) {
        error_msg("one or more input parameters missing.\n");
        return;
    }
	
    if(nrhs > 3) {
        error_msg("too many input parameters.\n");
        return;
    }
	
    const mxArray *X = prhs[0];    
    if(!mxIsDouble(X) | mxIsSparse(X) | mxIsComplex(X)) {
        error_msg("X must be a real full matrix of doubles.\n");
        return;
    }

    int X_M=mxGetM(X); // Each data point has M coordinates
	int X_N=mxGetN(X); // N data points
	double *X_pr=mxGetPr(X);

	
	int Y_M; // Each data point has M coordinates
	int Y_N; // N data points
	double *Y_pr;

	if (nrhs==3) {
		const mxArray *Y = prhs[1];    
		if(!mxIsDouble(Y) | mxIsSparse(Y) | mxIsComplex(Y)) {
			error_msg("Y must be a real full matrix of doubles.\n");
			return;
		}
		
		Y_M=mxGetM(Y); 
		Y_N=mxGetN(Y); 
		Y_pr=mxGetPr(Y);
		
		
		if (X_M!=Y_M) {
			error_msg("X and Y must have the same number of coordinates (same number of rows)\n");
		}
	}
	else {
		Y_M=X_M;
		Y_N=X_N;
		Y_pr=X_pr;
	}

	int NN;
	
	if ((mxGetM(prhs[nrhs-1])!=1) | (mxGetN(prhs[nrhs-1])!=1)){
		error_msg("NN must be an integer scalar\n");
		return;
	}

	GetInteger(prhs[nrhs-1],&NN);
	
	if (NN>X_N){
		char *msg[200];
		sprintf((char*)msg,"%d nearest neighbors requested. Only %d data points available. Using NN=%d nearest neighbors.",NN,X_N,X_N);
		mexWarnMsgTxt((char*)msg);
		NN=X_N;
	}
	
#ifdef DEBUG
	printf("Input parameters parsed successfully...\n");

	printf("X: M=%d  N=%d\n",X_M,X_N);
	printf("Y: M=%d  N=%d\n",Y_M,Y_N);
	printf("NN=%d\n",NN);
#endif
	
	mxArray *idx       = mxCreateDoubleMatrix(NN, Y_N, mxREAL);
	mxArray *distances = mxCreateDoubleMatrix(NN, Y_N, mxREAL);
	
    if ((idx==NULL) | (distances==NULL)) {
        error_msg("Failed to allocate output arrays...aborting\n");
        return;
    }
    
    double *idx_pr = mxGetPr(idx);
	double *distances_pr = mxGetPr(distances);

	double *row_dists=(double*)mxMalloc(X_N*sizeof(double));
	int    *row_idx=(int*)mxMalloc(X_N*sizeof(int));
    
    if ((row_dists==NULL) | (row_idx==NULL)) {
        error_msg("Failed to allocate auxiliry arrays...aborting\n");
        return;
    }
    
	
	for (int i=0; i<Y_N; i++) {				
		// compute all distances for point i		
		for (int j=0; j<X_N; j++) {
			
			register double sum = 0.0;
			register double temp;
			
			for(int k=0; k < X_M; k++) {
				temp = Y_pr[i*Y_M+k]-X_pr[j*X_M+k];
				sum += temp*temp;
			}
			
			row_idx[j]=j;
			row_dists[j]=sum;			
		}
		// sort the distances of row i and copy the largest NN  
		// distances and indices to the output array
		
		QuickSort(row_dists, row_idx, 0, X_N-1);
                
		for (int k=0; k<NN;k++) {
			idx_pr[i*NN+k]=(double)row_idx[k];
			distances_pr[i*NN+k]=row_dists[k];			
		}
		
	}
	
	
	// Sort each column of distances and idx.
	
	mxFree(row_idx);
	mxFree(row_dists);
	
	plhs[0] = idx;
	if(nlhs > 1)
		plhs[1] = distances;
	else
		mxDestroyArray(distances);		
}