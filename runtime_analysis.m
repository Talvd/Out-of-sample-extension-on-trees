%
% script for runtime analysis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

warning('off','all')
%runtime = (log(n))^2
k=5;
NN = 5;

times = [];

for N =100:100:2000
    %create the data for 2D samples

    angle = linspace(-pi,pi,N);
    points = [cos(angle);sin(angle)];

    delta = abs((angle(1)-angle(2))*.5);
    new_angle = linspace(-pi+delta,pi+delta,2*N-1);
    new_points =[cos(new_angle);sin(new_angle)];
    new_points = new_points(:,1:100);

    f =  cos(2*points);
    f = sqrt(sum(f.*f))';

    
    folder_name = ['runtime_analysis/data_2D_LMW_',num2str(N)];
    file_name = [folder_name,'/2D_data.mat'];

    done=0;    
    try
        load(file_name);
        done=1;
    catch err
        done=0;
    end

    if ~done
        T_data_partition = binary_tree(points, @binary_spectral_partition_med, k, @BruteForceNN, NN,0);

        [B, T] = Laplacian_multiwavelets_basis( points, k, T_data_partition, @BruteForceNN, NN , 0);
   
        mkdir(folder_name)
        cd(folder_name)
        save('2D_data', 'B', 'T');
        cd '../../'
    end
    tic  
    
    alpha = B'*f;
    f_new=zeros(100,1);

    for i=1:100
        p=new_points(:,i);
        [ext_B, T_ext] = LMW_extrapolation(points,p, T,B, NN);
        f_new(i) = ext_B(end,:)*alpha; 
    end

    t = toc;
    times = [times t];
   
end

nn= 100:100:2000;
figure
plot(nn,times);
hold
tt = log(nn);
c = tt(end)/times(end);
tt1 = tt/c;
plot(nn,tt1,'r');

tt2 = log(nn).^2;
c = tt2(end)/times(end);
tt1 = tt2/c;
plot(nn,tt1,'g');

c = nn(end)/times(end);
tt1 = nn/c;
plot(nn,tt1,'k');


