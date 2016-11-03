function [ V, E, Epsi] = GLE( points, NN_func, NN, k, func_eps, metric, percent)
%---------------------------------------------------------------
% Sparse Graph Laplacian Eigenvectors (GLE) implementation
%
% Input
%  points - the data points: d X N where d is the dimension
%  k - the number of required eigenvectors
%  NN - Number of nearest neighbours to consider in the affinity
%       construction matrix
%
% Output
%  V - N X k matrix. The eigenvectors are the coolumns.
%
%---------------------------------------------------------------
% Nir Sharon 29-10-12 and June 2013
%---------------------------------------------------------------

if length(points)>30000
    NN_func = @RANN;
    opts.maxit  = 4000;
else
    NN_func = @KNN_Matlab;
    opts.maxit  = 1000;
end
    

if nargin<7
    percent = [];
end

if nargin<6
    metric = [];
end
    
if (nargin<5)||isempty(func_eps)
    func_eps = @Set_Epsilon_V2;
end

if nargin<4
    k=6;
end

if (NN<=0)||(mod(NN,1)~=0)
	warning('Bad nearest neighbour parameter. NN is setted to defualt 5');
    NN = 5;
end

if size(points,2)<k
    warning('Too many eigenvectors or not enough data points');
    k = size(points,2);
end

if size(points,2)<(NN+1)
    NN = size(points,2)-1;
end

if size(points,2)==1
    V = points;
    return
end

% for 400X400 run
if length(points)>70
    func_eps = @Set_Epsilon_V3;
else
   func_eps = @Set_Epsilon_V2;
end

% I try to overcome the bug of repeated identical points (the D matrix must
% have non-zeros
MIN_Neighbours = 5;
if size(points,2)<MIN_Neighbours
    NN=size(points,2);
end

% if we choose NN==1 there are no neighbours (the case of single data point
% was already treated).
if NN==1
   NN=2;
end
   
[ind, dist] = NN_func(points,NN);

% if we given different metric
if ~isempty(metric)
    for i=2:(size(ind,1))
        for j=1:(size(ind,2))
            dist(i,j) = metric(points(:,j),points(:,ind(i,j)));
        end
    end
end

%%========= DEBUGGING PART ===================
% if (isequal(dist,zeros(size(dist))))
%     flag = 0;
% end
%%========= DEBUGGING PART ===================

% Epsi = median(dist(:));

%%======== Nir Code ============
Epsi = func_eps(ind, dist, 0, percent);  % auto select epsilon for the Graph Laplacian

[W,D]=Graph_Matrix(ind, dist, Epsi, 0);
%%======== Nir Code ============

% [W,D]=Graph_Matrix_Var(ind, dist, size(points,1));

D_minus_half = diag(sqrt(D.^(-1)));
S = D_minus_half*W;
S = S*D_minus_half;
S = .5*(S+S');

D_half = diag(sqrt(D));
% largest eigenvalues + corresponding eigenvectors
if k<=(size(points,2)/2)
    try
        % nX1 starting vector
%         opts.v0     = D_half*ones(numel(D),1);
        opts.isreal = 1;
        opts.issym  = 1;
        [V,E] = eigs(S,k,'lm',opts);
    catch  err % TO FIX !
        if (strcmp(err.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14'))
            opts.tol = 1e-3;
            [V,E]  = eigs(S, k, 'lr', opts);
            %[fV,E] = eig(full(S));
            [~, ind_e] = sort(diag(E),'descend');
            if k>size(points,2)
                k=size(points,2);
            end
            V = V(:,ind_e(1:k));
        else
            rethrow(err);
        end
        
    end
    
    [~, ind_e] = sort(diag(E),'descend');
    V = V(:,ind_e(1:k));
else
    [fV,E] = eig(full(S));
    [~, ind_e] = sort(diag(E),'descend');
    if k>size(points,2)
        k=size(points,2);
    end
    V = fV(:,ind_e(1:k));
end
V = D_minus_half*V;

if any(imag(V(:)))
    
    warning(['Complex GLE, eps is', num2str(Epsi),'length ', num2str(length(points))]);
    
    problem =1;
end

end

