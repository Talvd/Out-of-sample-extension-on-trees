/* Andrei Osipov (Yale University) */
#include <mex.h>
#include <stdio.h>
#include <math.h>

#undef VERBOSE
/* #define VERBOSE*/

void print_inner(int ion, char *str, double *a, int *ia, long int* la, int *n, int itype);
extern int tree_test__(long int *n, int *m, double *a,
		       int *numit, int *isuper, int *istat,
		       long int *near_convicts, double* dists_convicts,
		       int *k, double *w, long int *lenw);
extern int get_memory_size__(long int *n, int *m, int *k, int *numit, long int *nwords);

/*

syntax: rann(
   a[m][n], int32 inn[k][n], dnn[k][n], int32 numit, stats[100]
)

 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
  int m, n, m_inn, n_inn, k, m_dnn, n_dnn, i, numit, lenw, m_stats, n_stats;
  int nwords, isuper, istat;
  double *a, *pinn, *dnn, *w, *dnumit, *stats, *disuper, *distat;
  int *inn, *pnumit;
  */
  int numit, i, k_short, m_short, isuper, istat, nthreads, ion;
  long int m, n, m_inn, n_inn, k, m_dnn, n_dnn, m_stats, n_stats, nwords, lenw;
  int *pnumit, *pisuper, *pistat;
  double *dnumit, *stats, *a, *pinn, *dnn, *w, *disuper, *distat;
  long int *inn;
  int itimes[15];

  #ifdef VERBOSE
  mexPrintf("sizeof(long int) = %d\n", sizeof(long int));
  
  mexPrintf("nlhs = %d\n", nlhs);
  mexPrintf("nrhs = %d\n", nrhs);
  
  #endif


  if (nrhs != 7) {
    mexErrMsgTxt("usage: rann(a[m][n], int64 inn[k][n], dnn[k][n], int32 numit, int32 isuper, int32 istat, stats[100])");
  }

  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  
  #ifdef VERBOSE
  mexPrintf("m = %ld\n", m);
  mexPrintf("n = %ld\n", n);
  #endif

  m_inn = mxGetM(prhs[1]);
  n_inn = mxGetN(prhs[1]);

  if (n_inn != n) {
    mexErrMsgTxt("usage: rann(a[m][n], int64 inn[k][n], dnn[k][n], int32 numit, int32 isuper, int32 istat, stats[100])");
  }
  k = m_inn;
  
#ifdef VERBOSE
  mexPrintf("k = %d\n", k);
#endif

  m_dnn = mxGetM(prhs[2]);
  n_dnn = mxGetN(prhs[2]);

  if (m_dnn != k) {
    mexErrMsgTxt("usage: rann(a[m][n], int64 inn[k][n], dnn[k][n], int32 numit, int32 isuper, int32 istat, stats[100])");
  }
  if (n_dnn != n) {
    mexErrMsgTxt("usage: rann(a[m][n], int64 inn[k][n], dnn[k][n], int32 numit, int32 isuper, int32 istat, stats[100])");
  }

  m_stats = mxGetM(prhs[4]);
  n_stats = mxGetN(prhs[4]);
  if ( (m_stats != 1) || (n_stats != 100) ||
       mxIsSparse(prhs[4]) || mxIsComplex(prhs[4]) || 
       (!mxIsClass(prhs[4], "double")) ) {
    mexErrMsgTxt("fifth argument should be a real dense 1 x 100 matrix");
  }

  if (mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) || 
      (!mxIsClass(prhs[0], "double"))) {
	mexErrMsgTxt("first argument should be a real dense matrix");
      }

  if (mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) || 
      (!mxIsClass(prhs[1], "int64"))) {
	mexErrMsgTxt("second argument should be an int64 dense matrix");
      }

  if (mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]) || 
      (!mxIsClass(prhs[2], "double"))) {
	mexErrMsgTxt("third argument should be a real dense matrix");
      }

  dnumit = mxGetPr(prhs[3]);
  pnumit = (int *)dnumit;
  numit = *pnumit;
  
  #ifdef VERBOSE  
  mexPrintf("numit = %d\n", numit);
  #endif

  disuper = mxGetPr(prhs[5]);
  pisuper = (int *)disuper;
  isuper = *pisuper;
  
  #ifdef VERBOSE
  mexPrintf("isuper = %d\n", isuper);
  #endif

  distat = mxGetPr(prhs[6]);
  pistat = (int *)distat;
  istat = *pistat;
  
  #ifdef VERBOSE
  mexPrintf("istat = %d\n", istat);
  #endif

  stats = mxGetPr(prhs[4]);

  #ifdef VERBOSE
  mexPrintf("starting naive nearest neighbors\n");
  #endif

  a = mxGetPr(prhs[0]);
  pinn = mxGetPr(prhs[1]);
  inn = (long int *)pinn;
  dnn = mxGetPr(prhs[2]);
  
  #ifdef VERBOSE
  mexPrintf("a = %f, %f, %f, %f\n", *a, *(a+1), *(a+2), *(a+3));

  mexPrintf("inn = %ld, %ld, %ld, %ld\n", *inn, *(inn+1), *(inn+2), *(inn+3));

  mexPrintf("dnn = %f, %f, %f, %f\n", *dnn, *(dnn+1), *(dnn+2), *(dnn+3));

  mexPrintf("a = %f, %f, %f, %f\n", *a, *(a+1), *(a+2), *(a+3));
  #endif

  i = 4;
  /*  prin22_(a, &i); */
  
  #ifdef VERBOSE
  mexPrintf("after prin2\n");
  #endif
  
  nwords = 0;
  m_short = m;
  k_short = k;
  
  #ifdef VERBOSE
  mexPrintf("m_short = %d\n", m_short);
  mexPrintf("k_short = %d\n", k_short);
  #endif

  ion = 0;
  print_on_off__(&ion);

  nthreads = 4;
  get_threads_number_(&nthreads);
  
  #ifdef VERBOSE
  mexPrintf("nthreads = %d\n", nthreads);
  #endif
  get_memory_size_(&n, &m_short, &k_short, &numit, &isuper, &nthreads, &nwords);

  #ifdef VERBOSE
  mexPrintf("RANN needs to allocate array of %ld doubles\n", nwords);
  #endif

  lenw = nwords;
  w = (double *)mxMalloc(sizeof(double)*lenw);
  if (w == 0) {
    mexPrintf("Need memory for %ld real numbers\n", lenw);
    mexErrMsgTxt("Cannot allocate enough memory");
  }


  tree_test_(&n,&m_short,a,&numit,&isuper,&istat,inn,dnn,&k_short,w,&lenw);
  /*
    stats[0] = average square of distance to true nearest neighbors
    stats[1] = average square of distance to suspects
    stats[2] = average proportion of true nn among suspects
    stats[3] = memory required by RANN (in double words)
    stats[4] = time by iterations + supercharging
    stats[5] = time by statistics

    stats[10] = time of all_iter (sec)
    stats[11] = time of get_stats (sec)
    stats[12] = total of all_iter -> all_iter0
    stats[13] = total of all_iter -> all_iter0 -> one_iter
    stats[14] = total of all_iter -> all_iter0 -> one_iter -> centralize
    stats[15] = total of all_iter -> all_iter0 -> one_iter -> rotate
    stats[16] = total of all_iter -> all_iter0 -> one_iter -> one_iter0
    stats[17] = total of all_iter -> all_iter0 -> one_iter -> one_iter0 ->strbld
    stats[18] = total of all_iter -> all_iter0 -> one_iter -> one_iter0 ->boxes
    stats[19] = total of all_iter -> all_iter0 -> merge
    stats[20] = total of second_search
    stats[21] = total of second_search -> copy distances + indices
   */
  stats[0] = w[19];
  stats[1] = w[20];
  stats[2] = w[21];
  stats[3] = 1e0 * lenw;
  stats[4] = w[22];
  stats[5] = w[23];

  i = 24;
  prin2_("w = *", w, &i);  

  for (i = 0; i < 12; i++) {
    itimes[i] = w[29+i];
    stats[10+i] = w[29+i];
  }
  
  #ifdef VERBOSE
  i = 12;
  
  prinf_("itimes = *", itimes, &i);
  mexPrintf("TIMING (in seconds)\n"); 

  i = 1;
  prinf_("TIME for all_iter=*",itimes,&i);
  prinf_("TIME for get_stats=*",itimes+1,&i);

  prinf_("Total all_iter0=*",itimes+2,&i);
  prinf_("   -> one_iter=*", itimes+3,&i);
  prinf_("   -> -> centralize=*",itimes+4,&i);
  prinf_("   -> -> rotate=*",itimes+5,&i);
  prinf_("   -> -> one_iter0=*",itimes+6,&i);
  prinf_("   -> -> -> struct_bld=*",itimes+7,&i);
  prinf_("   -> -> -> boxes=*",itimes+8,&i);
  prinf_("   -> merge=*",itimes+9,&i);
  prinf_("Total second_search=*",itimes+10,&i);
  prinf_("   -> copy to dst,susp=*",itimes+11,&i);
  #endif

  mxFree(w);

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
