/* Andrei Osipov (Yale University) */
#include <mex.h>
#include <stdio.h>
#include <math.h>

#include <sys/stat.h>

/* Yoel suggested to add prototypes */
void prin22_(double* a, int* n);
void print_inner(int ion, char *str, double *a, int *ia, long int *la, int *n, int itype);
void prin2_(char *str, double *a, int *n);

extern int load_distances64__(int *k, long int *n, double *dsts);
extern int load_indices64__(int *k, long int *n, long int* idxs);

/*

syntax: rann_load_chunk(int64 idxs[k][n], dsts[k][n]), n is int64

*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int m, m_inn, k, m_dnn, i, numit, lenw, m_stats, n_stats;
  int nwords, isuper, istat;
  double *a, *pinn, *dnn, *w, *dnumit, *stats, *disuper, *distat;
  int *pnumit, *pisuper, *pistat;
  long int n_inn, n_dnn, n;
  long int *inn;



  mexPrintf("nlhs = %d\n", nlhs);
  mexPrintf("nrhs = %d\n", nrhs);

  if (nrhs != 2) {
    mexErrMsgTxt("usage: rann_load_chunk(int64 idxs[k][n], dsts[k][n])");
  }

  m_inn = mxGetM(prhs[0]);
  n_inn = mxGetN(prhs[0]);
  mexPrintf("m_inn = %d\n", m_inn);
  mexPrintf("n_inn = %ld\n", n_inn);

  m_dnn = mxGetM(prhs[1]);
  n_dnn = mxGetN(prhs[1]);
  mexPrintf("m_dnn = %d\n", m_dnn);
  mexPrintf("n_dnn = %ld\n", n_dnn);

  n = n_inn;
  k = m_inn;
  mexPrintf("k = %d\n", k);
 

  if (n_dnn != n) {
    mexErrMsgTxt("usage: rann_load_chunk(int64 idxs[k][n], dsts[k][n])");
  }
  if (m_dnn != k) {
    mexErrMsgTxt("usage: rann_load_chunk(int64 idxs[k][n], dsts[k][n])");
  }

  if (mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) || 
      (!mxIsClass(prhs[0], "int64"))) {
	mexErrMsgTxt("first argument should be an int32 dense matrix");
      }

  if (mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) || 
      (!mxIsClass(prhs[1], "double"))) {
	mexErrMsgTxt("third argument should be a real dense matrix");
      }

  pinn = mxGetPr(prhs[0]);
  inn = (long int *)pinn;
  dnn = mxGetPr(prhs[1]);

  check_indices_files(k,n);
  mexPrintf("starting loading indices from rann_idxs\n");
  load_indices64__(&k,&n,inn);
  mexPrintf("Done loading indices.\n");

  check_distances_files(k,n);
  mexPrintf("starting loading distances from rann_dsts\n");
  load_distances64__(&k,&n,dnn);
  mexPrintf("Done loading distances.\n");

}


void check_distances_files(int k, long int n) {
  return;
  struct stat st;
  char file_name[12];
  int i, j;
  long int ichunk, lens[10000], lentot;

  ichunk = 125000000;
  lentot = 0;
  for (i = 0; i < 9999; i++) {
    sprintf(file_name, "rann_dst%04d",i);
    mexPrintf("%s\n",file_name);
    
    if (stat(file_name,&st) != 0) goto lbl1;
    lens[i] = st.st_size;
    lentot = lentot + lens[i];

  }

  lbl1:
  mexPrintf("i = %d\n", i);
  mexPrintf("lentot = %ld\n", lentot);

  if (lentot != (8*n*k)) {
    mexPrintf("Total length should be %ld\n", 8*n*k);
    mexErrMsgTxt("Error: rann_dst have incorrect size.");
    return;
  }

  for (j = 0; j < (i-1); j++) {
    if (lens[j] != (8*ichunk)) {
      sprintf(file_name,"rann_dst%04d",j);
      mexPrintf("File %s has length %ld instead of %ld\n", file_name, lens[j], 8*ichunk);
      mexErrMsgTxt("Error: rann_dst have incorrect size.");
      return;
    }
  }
}

void check_indices_files(int k, long int n) {
  return;
  struct stat st;
  char file_name[12];
  int i, j;
  long int ichunk, lens[10000], lentot;

  ichunk = 125000000;
  lentot = 0;
  for (i = 0; i < 9999; i++) {
    sprintf(file_name, "rann_idx%04d",i);
    mexPrintf("%s\n",file_name);
    
    if (stat(file_name,&st) != 0) goto lbl1;
    lens[i] = st.st_size;
    lentot = lentot + lens[i];

  }

  lbl1:
  mexPrintf("i = %d\n", i);
  mexPrintf("lentot = %ld\n", lentot);

  if (lentot != (8*n*k)) {
    mexPrintf("Total length should be %ld\n", 8*n*k);
    mexErrMsgTxt("Error: rann_idx have incorrect size.");
    return;
  }

  for (j = 0; j < (i-1); j++) {
    if (lens[j] != (8*ichunk)) {
      sprintf(file_name,"rann_idx%04d",j);
      mexPrintf("File %s has length %ld instead of %ld\n", file_name, lens[j], 8*ichunk);
      mexErrMsgTxt("Error: rann_idx have incorrect size.");
      return;
    }
  }
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
