function [f_new] = run_nystrom_extension(f, points, new_points, eps)
%
% find the Unnormalized Graph Laplacian Eigenvectors and extend them to new points by nystrom extension
% to approximate the value of empirical function for the new points
%
% f - empirical function values
% points - the points to create the basis
% new_points - points to extend the basis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

if (nargin<4)||isempty(eps)
    eps = [];
end

NN = size(points,2);

[B,E, eps] = GLE_unnormalized( points, NN, size(points,2), eps);


alpha = B'*f;

f_new = [];
for p=new_points    
    ext_B = nystrom_extension(points, eps, p, B, E, 0);
    f_new =[f_new ext_B(end,:)*alpha];   
end

f_new = f_new';