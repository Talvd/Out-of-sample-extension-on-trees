% Andrei Osipov (Yale University)
function test_rann_save_chunk(a)

[m, n] = size(a)

unix('rm rann_save_chunk_fort.o');
unix('gfortran -c -fPIC -fsecond-underscore rann_save_chunk_fort.F');
mex rann_save_chunk.c rann_save_chunk_fort.o -lgfortran -largeArrayDims
rann_save_chunk(a);
