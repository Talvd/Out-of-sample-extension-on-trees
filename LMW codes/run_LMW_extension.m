function [f_new] = run_LMW_extension(f, points, new_points, k,d, extension_func)
%
% create LMW basis and extend this basis to new points to approximate
% the value of empirical function for the new points
%
% f - empirical function values
% points - the points to create the basis
% new_points - points to extend the basis
% k - the LMW parameter
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

if (nargin<6)||isempty(extension_func)
    extension_func = @quaddric_fit_extension;
end

NN = 5;

%create the LMW basis
%--------------------------------------------------------------------------

folder_name = ['data_' num2str(d) 'D_LMW_',num2str(k)];
file_name = [folder_name,'/' num2str(d) 'D_data.mat'];


try
    load(file_name);
    done=1;
catch 
    done=0;
end

if ~done
    % binary tree based on spectral partition
    T_data_partition = binary_tree(points, @binary_spectral_partition_med, k, @BruteForceNN, NN,0);

    [ B, T ] = Laplacian_multiwavelets_basis( points, k, T_data_partition, @BruteForceNN, NN , 0);
   
    mkdir(folder_name)
    cd(folder_name)
    save([num2str(d) 'D_data'], 'B', 'T');
    cd '../'
end

%--------------------------------------------------------------------------

alpha = B'*f;
f_new=[];

for p=new_points    
    [ext_B, T_ext] = LMW_extrapolation(points,p, T,B, NN, extension_func);
    f_new = [f_new ext_B(end,:)*alpha]; 
end

f_new = f_new';
