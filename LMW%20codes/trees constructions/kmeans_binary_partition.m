function [ part1, part2 ] = kmeans_binary_partition( points, ~, ~ )
%
% partition based on Matlab's kmean implementation
%
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

if n==0
    error('empty array to subdivide');
end

ind_arr = 1:n;

if n>1
    partition = kmeans(points',2);
else
    partition = 1;
end

part1 = ind_arr(partition==1);
part2 = ind_arr(partition==2);


end

