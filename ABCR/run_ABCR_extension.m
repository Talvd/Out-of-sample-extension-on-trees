function [f_new] = run_ABCR_extension(f, points, new_points)
%
% create ABCR basis and extend this basis to new points to approximate
% the value of empirical function for the new points
%
% f - empirical function values
% points - the points to create the basis
% new_points - points to extend the basis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

l = 4;
k = size(points,2)/(2^l);

%create the ABCR basis
%--------------------------------------------------------------------------

folder_name = ['data_1D_ABCR_',num2str(k)];
file_name = [folder_name,'/1D_data.mat'];


try
    load(file_name);
    done=1;
catch 
    done=0;
end

if ~done
    B=ABCR_Version4_Overall( points,k,l);

    mkdir(folder_name)
    cd(folder_name)
    save('1D_data', 'B');
    cd '../'
end

%--------------------------------------------------------------------------

alpha = B'\f;
f_new=[];

for p=new_points    
    [ext_B] = ABCR_extrapolation(B, points, p, k);
    f_new = [f_new ext_B(:,end)'*alpha]; 
end

f_new = f_new';