function [ B, T ] = Laplacian_multiwavelets_basis( points, k, T, NN_func, NN , plotting, func_eps, metric)
%---------------------------------------------------------------
% Construct a binary hierarchy tree using binary partition method (binary_partition)
% We use the tree structure
%
% Input
%  points - the data points: d X N where d is the dimension
%  k - the number of minimal points per partition subset.
%  NN_func - the nearest neighbours function in use.
%  NN - Number of nearest neighbours to consider in the affinity
%       construction matrix
%  plotting - whether to plot or not the trees along the process
%
% Output
%  B - the LMW basis
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

if nargin<4
    error('Not enough arguments')
end

if nargin<5
    NN = 8;
end

if nargin<6
    plotting = 0;
end

if nargin<7
    func_eps = [];
end

if nargin<8
    metric = [];
end

if plotting
    figure(1);
    plot1(T,@node2str_nodenum);
    title('Partition: number of elements');
end

% first step
tic
T = Approximation_Spaces_V2( points, k, T, NN_func, NN, func_eps, metric  );
t1 = toc;
if plotting
    figure(2);
    plot1(T,@node_rank);
	title('Local rank, nested spaces');
end

% second step
tic
T = Fast_Orthogonalization(T);
t2 = toc

if plotting
    figure(3);
    plot1(T,@node_rank_second);
	title('Rank of complement spaces, except leafs');
end

% packing the basis into a sparse matrix
B = tree2basis(T);

if plotting
    figure(4);
    imagesc(B);
    m = OrthoCheck(B');
	title(['the resulted basis, orthogonality of ',num2str(m)]);
end

end

