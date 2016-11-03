function [approx_f ext_f S] = singleScale(f,in_dists,out_dists,ep)
% INPUT
% f        : an nx1 vector of the function values we interpulate 
%            that is f0[i] = f(x_i) for x_i in X=(x_1,...,x_n)
% in_dists : an nxn symmetric matrix describing the distances between the
%            points f0 is defined on ; in_dist(i,j) = dist(x_i,x_j)
% out_dists: an mxn matrix containing the distances for the m new points
%            which we wish to approximate. out_dist(i,j) = dist(y_i,x_j) for 1=<i=<m
% ep       : the width of the gaussians we use in the approximation

% OUTPUT
% approx_f :  an n x 1 vector containing the f values for the approximated f
% S        :  an n x 1 vector zero/one matrix containing the sampled
%              x_i : S(i) is one if x_i was sampled 
% ext_f    : Extanded F. m x 1 vector containing the approximated values of y_i

EIGVAL_TRESH = 0.1; % how many eigenvalues to use

G = exp(-in_dists.^2/ep);
k = rank(G/norm(G),EIGVAL_TRESH);
S = deterministicInterpolativeDecomposition(G,k);
coefs = G(:,S)\f;   %  i.e. coefs = inverse(G(:,S)) * f
approx_f = G(:,S)*coefs;
ext_f = exp(-(out_dists(:,S)).^2/ep)*coefs; 
