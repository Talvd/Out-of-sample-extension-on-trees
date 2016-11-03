N = 150;
d = 2;
k = 6;
NN = 7;

points = rand(d,N);
T = binary_tree(points, @binary_spectral_partition, k, NN);