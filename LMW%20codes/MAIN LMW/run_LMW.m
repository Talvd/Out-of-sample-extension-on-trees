%---------------------------------------------------------------
% Run the main LMW function "Laplacian_multiwavelets_basis"
% 
% k - the LMW parameter
% N - number of points to create
% NN - nearest neighbour parameter
% d - data dimension
%                                                                       
% Other Variables
% 
% B,B_haar,B_lap - the bases, N X N matrices. Elements are given in columns
% f - the function. Given in column, N X 1.
% points - data points.
% 
%---------------------------------------------------------------
% Nir Sharon, June 2013.
%---------------------------------------------------------------

%----------------- PARAMETER FOR THE RUN -----------------------
N = 200;
d = 3;
k = 5;
NN = 7;
binarytree = 1;
%------------------ CREATING THE DATA --------------------------
points = rand(d,N);


points = 2*points-ones(d,N);
norms  = sqrt(sum(points.*points));
points = points*diag(norms.^(-1));

scatter3(points(1,:),points(2,:),points(3,:))

% space = 1/N;
% points = 0:space:1;

% % N=120;
% % k=1;  % HAAR
% % points = linspace(-1,1,N);

%points = sort(rand(1,N));

%-------------- CONSTRUCTING THE TREE --------------------------
if binarytree
    % binary tree based on spectral partition
    T_data_partition = binary_tree(points, @binary_spectral_partition_med, k, @BruteForceNN, NN);
else
    % k-means partition
    max_degree_of_tree_node = 4;
    T_data_partition = hierarchy_tree(points, @kmeans_based_partition, k, max_degree_of_tree_node);
end
%---------------------- THE RUN ---------=====------------------
[B, T] = Laplacian_multiwavelets_basis( points, k, T_data_partition, @BruteForceNN, NN , 1);

norm(full(B))-1



