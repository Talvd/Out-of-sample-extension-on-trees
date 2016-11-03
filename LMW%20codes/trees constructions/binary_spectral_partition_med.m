function [ part1, part2 ] = binary_spectral_partition_med( points, NN_func, NN )
% partition by the sign of the second eigenvector of
% the graph Laplacian.
% Input
%  points - the data points: d X N where d is the dimension
%  NN - Number of nearest neighbours to consider in the affinity
%       construction matrix
%
% Output
%  part1,part2 - two sets of indices of the partition of "points"
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

n = size(points,2);
ind_arr = 1:n;


if n>1
     %v = GLE(points,NN_func,NN,2); % using @Set_Epsilon_V2
    v = GLE(points,NN_func,NN,2,@Set_Epsilon_V3);
    med = median(v(:,2));
    part1 = ind_arr(v(:,2)<med);
    part2 = ind_arr(v(:,2)>=med);
else
    part1 = 1;
    part2 = [];
end
end

