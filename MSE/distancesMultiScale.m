function [AF S EF errors] = distancesMultiScale(f0, in_dists, out_dists, err ,init_ep)
% Multi Scale Data Sampliing and Function Extension
% This function can estimate new values of a function f0 at given points based on distances
% from a know training set which we know the f0 values.
% The method can be viewed in 
% http://www.cs.tau.ac.il/~amir1/PS/MSE-ACHA-P.pdf

% INPUT
% f0        : an nx1 vector of the function values we interpulate 
%             that is f0[i] = f(x_i) for x_i in X=(x_1,...,x_n)
% in_dists  :  an nxn symmetric matrix describing the distances between the
% points f0 is defined on ; in_dist(i,j) = dist(x_i,x_j)
% out_dists : an mxn matrix containing the distances for the m new points
% which we wish to approximate. out_dist(i,j) = dist(y_i,x_j) for 1=<i=<m
% err       : the addmisable approximation error
% init_ep   : the initial width of the gaussians we use in the
% approximation

% OUTPUT
% AF        :  an n x iterations array containing the f values for each
% iteration
% S         :  an n x iterations zero/one matrix containing the sampled
%              x_i : S(i,j) is one if x_i was sampled in iteration j
% EF        : Extanded F. m x iteration matrix containing the approximated values of y_i
% errors    : 1 x iteration matrix containing the average L2 error

if ( (size(f0,2) ~=1)) || (length(f0) ~= length(in_dists)) || ...
         (size(in_dists,1) ~= size(in_dists,2))...
         || (size(out_dists,2) ~= size(in_dists,1))
     error('distancesMultiScale(): wrong size of input values.');
end

EPS_SHRINK_FACTOR = 4;
MAX_ITERATIONS = 100;

AF = [];
EF = [];
S = [];
EXT_f = 0;
APPROX_f = 0;
ep =init_ep;
f = f0;
Err = inf;
errors = [];
lf = length(f);
ls = 0;
count_iterations=0;

while ( Err > err && ls<lf && count_iterations<MAX_ITERATIONS)
    [approx_f ext_f s] = singleScale(f,in_dists,out_dists,ep);
    ls = length(s);
    v = zeros(length(f),1);
    v(s)=1;
    S = [S v];
    APPROX_f = APPROX_f + approx_f;
    AF = [AF APPROX_f];
    EXT_f = EXT_f + ext_f;
    EF = [EF EXT_f];    
    delta = f0- APPROX_f;
    f = delta;
    Err = norm(delta)/lf;
    errors = [errors Err];    
    ep = ep/EPS_SHRINK_FACTOR;
    count_iterations=count_iterations+1;        
end
