function [ V, E, Epsi] = GLE_unnormalized( points, NN, k , eps)
%---------------------------------------------------------------
% Unnormalized Sparse Graph Laplacian Eigenvectors (GLE) implementation
%
% Input
%  points - the data points: d X N where d is the dimension
%  k - the number of required eigenvectors
%  NN - Number of nearest neighbours to consider in the affinity
%       construction matrix
%
% Output
%  V - N X k matrix. The eigenvectors are the columns.
%

if length(points)>30000
    NN_func = @RANN;
    opts.maxit  = 4000;
else
    NN_func = @KNN_Matlab;
    opts.maxit  = 1000;
end
    
func_eps = @Set_Epsilon_V2;

if nargin<3
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

[ind, dist] = NN_func(points,NN);

if nargin<4 || isempty(eps)
    Epsi = func_eps(ind, dist, 0);  % auto select epsilon for the Graph Laplacian
else
    Epsi = eps;
end

[W,~]=Graph_Matrix(ind, dist, Epsi, 0);

S=W;
% largest eigenvalues + corresponding eigenvectors
if k<=(size(points,2)/2)
    try
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
    [fE, ind_e] = sort(diag(E),'descend');
    if k>size(points,2)
        k=size(points,2);
    end
    V = fV(:,ind_e(1:k));
    E = diag(fE);
end

if any(imag(V(:)))
    
    warning(['Complex GLE, eps is', num2str(Epsi),'length ', num2str(length(points))]);
    
    problem =1;
end

end

