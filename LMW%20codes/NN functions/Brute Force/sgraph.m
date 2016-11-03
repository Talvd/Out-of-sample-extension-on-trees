function [W,D]=sgraph(X,sigma,NN,threshold)
%
% Create a symmetric sparse graph where
%
%   W(i,j)=exp(-norm(X(i)-X(j))^2/sigma)
%
% The matrix W is symmetrized, but no normalization is applied.
% The matrix D contains the row sums of W for future normalization.
%
% sigma can be a vector, in which case, each point has a different sigma.
%
% If threshold is different from 0, edges that are smaller than the
% threshold are set to zero.
%
% NN is the number of nearest neighbors of each point.
%
% Default: sigma=1
%          NN=10
%
% Yoel Shkolnisky 31/08/2006

if nargin<4
    threshold=0;
end

if nargin<3
    NN=10;
    fprintf('Using default value NN=%d\n',NN);
end

if nargin<2
    sigma=1;
    fprintf('Using default value sigma=%d\n',sigma);
end

n=size(X,1);
atria=nn_prepare(X);

% The commented line should work but for some reason it is not, so use the
% following:
[index,distance]=nn_search(X,atria,[1:n].',NN,-1,0.0);

if prod(size(sigma))>2
    sigma=repmat(sigma,1,NN);
end

S=exp(-distance.^2./sigma);
S=S.';
S=S(:);

if threshold>0
    S(find(S<threshold))=0;
end

I=zeros(NN,n);
for k=1:NN
    I(k,:)=1:n;
end
I=I(:);

J=index.';
J=J(:);
W=sparse(I,J,S(:),n,n,n*NN);
W=(W+W.')/2;

D=sum(W);
