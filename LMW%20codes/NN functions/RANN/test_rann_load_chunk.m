% Andrei Osipov (Yale University)
function test_rann_load_chunk(n,k)


unix('rm rann_save_chunk_fort.o');
unix('gfortran -c -fPIC -fsecond-underscore rann_save_chunk_fort.F');
mex rann_load_chunk.c rann_save_chunk_fort.o -lgfortran -largeArrayDims

idxs = int64(zeros(k,n));
dsts = zeros(k,n);

rann_load_chunk(idxs,dsts);

idxs
dsts
