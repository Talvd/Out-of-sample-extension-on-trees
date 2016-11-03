% Andrei Osipov (Yale University)
function test_rann_shell_par(n,m,k,nit,isuper,istat)

file_diary = ['diary_rann_shell_par_n', int2str(n), '_m', int2str(m), ...
    '_k', int2str(k)];
delete(file_diary);
diary(file_diary);

diary on;

% delete files
unix('rm rann_pts0*');
unix('rm rann_idx0*');
unix('rm rann_dst0*');

% to be able to save and load points
unix('rm rann_save_chunk_fort.o');
unix('gfortran -c -fPIC -fsecond-underscore rann_save_chunk_fort.F');
mex rann_save_chunk.c rann_save_chunk_fort.o -lgfortran -largeArrayDims

% to be able to run the algorithm
unix('./par_shell');

% to be able to read indices and distances
unix('rm rann_save_chunk_fort.o');
unix('gfortran -c -fPIC -fsecond-underscore rann_save_chunk_fort.F');
mex rann_load_chunk.c rann_save_chunk_fort.o -lgfortran -largeArrayDims

randn('seed',0);
X = randn(m, n);
[m, n] = size(X)
m
n
nit

idxs = int64(zeros(k, n));
dsts = zeros(k, n);

% save points into rann_dir
rann_save_chunk(X)
clear('X');

% run algorithm
disp('Running rann_shell_par...');
told = now;
disp(['Started at ', datestr(now)]);

shell_str = ['./rann_shell_par ', int2str(n), ' ', int2str(m), ' ', int2str(k), ...
    ' ', int2str(nit), ' ', int2str(isuper), ' ', int2str(istat)]
ier = unix(shell_str)

disp(['Finished at ', datestr(now)]);
seconds = (now-told)*24*60*60

ier

if (ier == 0)
   
    % load indices from rann_idxs and distances from rann_dsts
    rann_load_chunk(idxs,dsts);
    
    disp('Start saving idx_par, dst_par...');
    file_idx = 'idx_par';
    file_dst = 'dst_par';
    idxs = double(round(idxs));
    save(file_idx, 'idxs', '-ascii');
    save(file_dst, 'dsts', '-ascii');
    disp('Done saving.');

    disp(['n=', int2str(n), ',m=',int2str(m),',k=',int2str(k),',nit=',int2str(nit)]);
    disp(['isuper=', int2str(isuper),',istat=',int2str(istat)]);
end

diary off;

end
