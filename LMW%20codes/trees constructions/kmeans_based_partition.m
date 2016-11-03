function [ p_cell ] = kmeans_based_partition( points, max_deg )
% partition by the sign of the second eigenvector of
% the graph Laplacian.
% Input
%  points - the data points: d X N where d is the dimension
%  max_deg - Number of maximal partitions
%       
% Output
%  part1,part2 - two sets of indices of the partition of "points"
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

n = length(points);
ind_arr = 1:n;

if nargin<2
    max_deg = 6 ;
    warning(['maximal degree is setted to ', int2str(max_deg)])
end

% choosing the clustering parameter
most = max_deg; 
msr = zeros(1,most);
msr(1) = inf;

for par=2:most
  ind = kmeans(points',par);
  s = silhouette(points',ind);
  msr(par) = mean(s);
end

[~, choose] = min(msr);

% the partition
partition = kmeans(points',choose);
p_cell = cell(choose,1);
for par=1:choose
    p_cell{par} = ind_arr(partition==par);
end

end

