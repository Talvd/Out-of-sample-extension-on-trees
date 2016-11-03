% Andrei Osipov (Yale University)
function test_points64(n,m,k,nit,isuper,istat)

file_diary = ['diary_test_points64_k', int2str(k)];
delete(file_diary);
diary(file_diary);

diary on;

unix('./rann64');
mex rann64.c rann_par_fort64.o corrand4.o random_transf4.o durdec.o peter_heapsort_rea64.o peter_heapsort_rea_part64.o peter_heapsort_inte64.o -lgfortran -largeArrayDims -lgomp

randn('seed',0);
X = randn(m, n);
[m, n] = size(X)
m
n
nit

idxs = int64(zeros(k, n));
dsts = zeros(k, n);

numit = int32(nit);
stats = zeros(1, 100);
iisuper = int32(isuper);
iistat = int32(istat);
disp('Running rann64...');
told = now;
disp(['Started at ', datestr(now)]);
rann64(X, idxs, dsts, numit, stats, iisuper, iistat);
disp(['Started at ', datestr(told)]);
disp(['Finished at ', datestr(now)]);
seconds = (now-told)*24*60*60


stats(1:3)
stats(4)
size(X)
size(idxs)
size(dsts)

id = idxs(:,1:10)
ds = dsts(:,1:10)

file_idx = ['idx_points64_k', int2str(k)]
file_dst = ['dst_points64_k', int2str(k)]
ids2 = double(round(idxs));
save(file_idx, 'ids2', '-ascii');
save(file_dst, 'dsts', '-ascii');

disp(['n=', int2str(n), ',m=',int2str(m),',k=',int2str(k),',numit=',int2str(numit)]);
disp(['isuper=', int2str(isuper),',istat=',int2str(istat)]);
disp(['average distance to true nn = ', num2str(stats(1),'%e')]);
disp(['average distance to suspects = ', num2str(stats(2),'%e')]);
disp(['ratio = ', num2str(stats(2)/stats(1),'%e')]);
disp(['ratio - 1 = ', num2str((stats(2)/stats(1)) - 1,'%e')]);
disp(['proportion of true nn among suspects = ', num2str(stats(3),'%e')]);
disp(['Memory required by the algorithm, in double words = ', num2str(stats(4),'%e')]);
disp(['Time by rann64, in seconds = ', num2str(seconds)]);
disp(['Time by iterations and supercharging only = ', num2str(stats(5))]);
disp(['Time by statistics = ', num2str(stats(6))]);

diary off;

end