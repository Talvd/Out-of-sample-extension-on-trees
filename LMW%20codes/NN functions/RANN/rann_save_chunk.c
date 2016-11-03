/* Andrei Osipov (Yale University) */
#include <mex.h>
#include <stdio.h>
#include <math.h>

void prin22_(double* a, int* n);
void print_inner(int ion, char *str, double *a, int *ia, long int *la, int *n, int itype);
void prin2_(char *str, double *a, int *n);

extern int save_points64__(int *m, long int *n, double *a);

/*

syntax: rann_save_chunk(a[m][n]), where n = int64

*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int m, m_inn, n_inn, k, m_dnn, n_dnn, i, numit, lenw, m_stats, n_stats;
  int nwords, isuper, istat;
  double *a, *pinn, *dnn, *w, *dnumit, *stats, *disuper, *distat;
  int *inn, *pnumit, *pisuper, *pistat;
  long int n;



  mexPrintf("nlhs = %d\n", nlhs);
  mexPrintf("nrhs = %d\n", nrhs);

  if (nrhs != 1) {
    mexErrMsgTxt("usage: rann_save_chunk(a[m][n])");
  }

  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  mexPrintf("m = %d\n", m);
  mexPrintf("n = %ld\n", n);

  if (mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) || 
      (!mxIsClass(prhs[0], "double"))) {
	mexErrMsgTxt("first argument should be a real dense matrix");
      }

  mexPrintf("starting saving points into a file\n");

  a = mxGetPr(prhs[0]);

  mexPrintf("a = %f, %f, %f, %f\n", *a, *(a+1), *(a+2), *(a+3));

  save_points64__(&m,&n,a);

  mexPrintf("points saved\n");
}


/*

  PRINTING FUNCTIONS TO BE CALLED FROM FORTRAN

 */

void prin22_(double* a, int* n) {

  /*  mexPrintf("a = \n"); */
  int i, ix;

  i = 0;
  ix = 0;
  for (i = 0; i<*n; i++) {
    mexPrintf("  %11.4le", *(a+i));
    ix = ix + 1;
    if (ix == 6) {
      ix = 0;
      mexPrintf("\n");
    }
  }
  if (ix > 0) {
    mexPrintf("\n");
  }

}


void cprina_(char *str) {

  int i;
  for (i = 0; i < 100; i++) {
    if (str[i] == '*') goto jump1;
    mexPrintf("%c", str[i]);
  }
 jump1:

  mexPrintf("\n");

}


void cprin2_(char *str, double *a, int* n) {
  cprina_(str);  
  prin22_(a, n);
}

void cprinf_(char *str, int *a, int *n) {
  int i, ix;

  cprina_(str);
  

  i = 0;
  ix = 0;
  for (i = 0; i<*n; i++) {
    mexPrintf(" %7d", *(a+i));
    ix = ix + 1;
    if (ix == 10) {
      ix = 0;
      mexPrintf("\n");
    }
  }
  if (ix > 0) {
    mexPrintf("\n");
  }


}

void cprinl_(char *str, long int *a, int *n) {
  int i, ix;

  cprina_(str);
  

  i = 0;
  ix = 0;
  for (i = 0; i<*n; i++) {
    mexPrintf(" %7ld", *(a+i));
    ix = ix + 1;
    if (ix == 10) {
      ix = 0;
      mexPrintf("\n");
    }
  }
  if (ix > 0) {
    mexPrintf("\n");
  }


}


/*

  SECOND GENERATION PRINTING FUNCTIONS

 */


void print_on_off__(int *ion) {
  char *str;
  double *a;
  int *ia;
  int *n;
  long int *la;
  print_inner(*ion, str, a, ia, la, n, 1);
}

void prina_(char* str) {
  double *a;
  int *ia;
  int *n;
  long int *la;
  print_inner(0, str, a, ia, la, n, 2);
}

void prin2_(char* str, double* a, int *n) {
  int *ia;
  long int *la;
  print_inner(0, str, a, ia, la, n, 3);
}

void prinf_(char* str, int* ia, int *n) {
  double *a;
  long int *la;
  print_inner(0, str, a, ia, la, n, 4);
}

void prinl_(char* str, long int* la, int *n) {
  double *a;
  int *ia;
  print_inner(0, str, a, ia, la, n, 5);
}

void print_inner(int ion, char *str, double *a, int *ia, long int *la, int *n, int itype) {

  static int iprint_on;

  /* PRINT ON / OFF */
  if (itype == 1) {
    iprint_on = ion;
    return;
  }

  if (iprint_on == 0) return;

  /* PRINT STRING */
  if (itype == 2) {
    cprina_(str);
  }

  /* PRINT DOUBLE */
  if (itype == 3) {
    cprin2_(str, a, n);
  }

  /* PRINT INT */
  if (itype == 4) {
    cprinf_(str, ia, n);
  }

  /* PRINT LONG INT */
  if (itype == 5) {
    cprinl_(str, la, n);
  }


}
