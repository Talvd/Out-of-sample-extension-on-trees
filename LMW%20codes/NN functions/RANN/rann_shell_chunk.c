/* Andrei Osipov (Yale University) */
#include <stdio.h>
#include <stdlib.h>

#include <sys/stat.h>

/*

  COMPILE by using rann_sh

  USAGE: rann_shell_chunk n m k numit isuper istat

  NEEDS: file rann_pts* 
     written by test_rann_save_chunk.m, using
                rann_save_chunk_fort.F, rann_save_chunk.c

  this file is of length 8*n*m bytes
  it contains the points to run the RANN on

  OUTPUT: two binary files
     rann_dst*
        size: 8*k*n bytes, contains distances to suspects
     rann_idx*
        size: 8*k*n bytes, contains indices of the suspects (int64)

  ERROR CODES: see below

  TODO: file rann_stat
        (probably a small ascii file to save statistics)

  HOW TO USE FROM MATLAB:
  1) Call test_rann_save_chunk(a), where a[m][n] are the points
  2) Call ier = unix('rann_shell_chunk 5 3 2 3 1 1'), where ier = error code
  3) TODO: Call test_rann_load_chunk(n,k)

 */



/*
  ERROR CODES

  1
     must be 6 arguments, n m k numit isuper istat
  2
     n is not a positive integer
  3
     m is not a positive integer
  4
     k is not a positive integer
  5
     numit is not a positive integer
  6
     isuper is not a positive integer
  7
     istat is not a positive integer
  8 
     k is greater than or equal to n
  9
     cannot read the file rann_dir64
  10
     the files rann_pts* must be of total length 8*n*m bytes
  11
     cannot allocate enough memory

 */

void main(int argc, char **argv)
{
  int i, m, k, numit, isuper, istat, nm_short, nk_short, ier, ion;
  long int n, lfile, nm, nk, nwords;
  double *a, *dsts, *w;
  long int *idxs;

  struct stat st;

  /*
  ion = 1;
  print_on_off__(&ion);
  */

  printf("argc = %d\n", argc);
  for (i = 0; i < argc; i++) {
    printf("arg[%d]=%s\n", i, argv[i]);
  }

  if (argc != 7) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    exit(1);
  }

  n = atol(argv[1]);
  if (n < 1) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("n must be a positive integer\n");
    exit(2);
  }
  printf("n = %7ld\n", n);

  m = atoi(argv[2]);
  if (m < 1) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("m must be a positive integer\n");
    exit(3);
  }
  printf("m = %d\n", m);

  k = atoi(argv[3]);
  if (k < 1) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("k must be a positive integer\n");
    exit(4);
  }
  printf("k = %d\n", k);

  numit = atoi(argv[4]);
  if (numit < 1) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("numit must be a positive integer\n");
    exit(5);
  }
  printf("numit = %d\n", numit);

  isuper = atoi(argv[5]);
  if (istat < 0) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("isuper must be a non-negative integer\n");
    exit(6);
  }
  printf("isuper = %d\n", isuper);

  istat = atoi(argv[6]);
  if (istat < 0) {
    printf("usage: rann_shell_chunk n m k numit isuper istat\n");
    printf("k must be a non-negative integer\n");
    exit(7);
  }
  printf("istat = %d\n", istat);

  if (k >= n) {
    printf("k cannot be greater than or equal to n\n");
    exit(8);
  }

  /* check the validity of rann_pts* files */
  check_point_files(m, n, &ier);
  if (ier > 0) {
    exit(ier);
  }


  printf("Allocating %7ld real words for the points...\n", n*m);
  a = (double *)malloc(sizeof(double)*m*n);
  if (a == 0) {
    printf("Cannot allocate enough memory\n");
    exit(11);
  }
  printf("Allocation complete.\n");

  printf("Reading points from rann_dir64...\n");
  load_points64__(&m,&n,a);

  nm = n*m;
  if (nm > 1000) nm = 1000;
  nm_short = nm;
  prin2_("a = *", a, &nm_short);

  printf("Allocating %7ld real words for the distances...\n", n*k);
  dsts = (double *)malloc(sizeof(double)*k*n);
  if (dsts == 0) {
    printf("Cannot allocate enough memory\n");
    exit(11);
  }
  printf("Allocation complete.\n");


  printf("Allocating %7ld integer words for the indices...\n", n*k);
  idxs = (long int *)malloc(sizeof(long int)*k*n);
  if (idxs == 0) {
    printf("Cannot allocate enough memory\n");
    exit(11);
  }
  printf("Allocation complete.\n");
  

  /* Determine how much memory is required */
  nwords = 0;
  get_memory_size__(&n, &m, &k, &numit, &nwords);
  printf("Allocating %7ld real words for the working array\n", nwords);

  w = (double *)malloc(sizeof(double)*nwords);
  if (w == 0) {
    printf("Cannot allocate enough memory\n");
    exit(11);
  }
  printf("Allocation complete.\n");

  printf("Running rann64...\n");
  tree_test__(&n,&m,a,&numit,&isuper,&istat,idxs,dsts,&k,w,&nwords);
  printf("Done running rann64.\n");

  nk = n*k;
  if (nk > 1000) nk = 1000;
  nk_short = nk;
  prin2_("dsts = *", dsts, &nk_short);
  prinl_("idxs = *", idxs, &nk_short);

  printf("Avg dist to true nn = %f\n", w[19]);
  printf("Avg dist to suspects = %f\n", w[20]);
  printf("Ratio = %f\n", w[20]/w[19]);
  printf("Ratio-1 = %f\n", w[20]/w[19]-1.0);
  printf("Proportion of true nn among suspects = %f\n", w[21]);
  printf("Memory required, in bytes = %7ld\n", 8*n*m+8*n*k+8*n*k+8*nwords);
  printf("Time by iterations and supercharging only = %f\n", w[22]);
  printf("Time by statistics = %f\n", w[23]);

  printf("Saving distances to rann_dst...\n");
  save_distances64__(&k,&n,dsts);
  printf("Done saving distances.\n");
  printf("Saving indices to rann_idx...\n");
  save_indices64__(&k,&n,idxs); 
  printf("Done saving indices.\n");

  /*
   now test reading
   */
  load_distances64__(&k,&n,dsts);
  prin2_("dsts = *", dsts, &nk_short);
  load_indices64__(&k,&n,idxs);
  prinl_("idxs = *", idxs, &nk_short);

  free(a);
  free(dsts);
  free(idxs);
  free(w); 

  exit(0);

}


void check_point_files(int m, long int n, int *ier) {
  struct stat st;
  char file_name[12];
  int i, j;
  long int ichunk, lens[10000], lentot;

  ichunk = 125000000;
  lentot = 0;
  for (i = 0; i < 9999; i++) {
    sprintf(file_name, "rann_pts%04d",i);
    printf("%s\n",file_name);
    
    if (stat(file_name,&st) != 0) goto lbl1;
    lens[i] = st.st_size;
    lentot = lentot + lens[i];

  }

  lbl1:
  printf("i = %d\n", i);
  printf("lentot = %ld\n", lentot);

  if (lentot != (8*n*m)) {
    printf("Total length should be %ld\n", 8*n*m);
    *ier = 19;
    return;
  }

  for (j = 0; j < (i-1); j++) {
    if (lens[j] != (8*ichunk)) {
      sprintf(file_name,"rann_pts%04d",j);
      printf("File %s has length %ld instead of %ld\n", file_name, lens[j], 8*ichunk);
      *ier = 20;
      return;
    }
  }
  *ier = 0;
}


