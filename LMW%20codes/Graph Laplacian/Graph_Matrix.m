function [W,D]=Graph_Matrix(index,distance,sigma,threshold)
%
% Create a symmetric sparse graph matrix where
%
%   W(i,j)=exp(-norm(X(i)-X(j))^2/sigma)
%
% Input
%
%  index,distance - are the info for the nearest neighbors of each point.
%
%  sigma can be a vector, in which case, each point has a different sigma.
%  default: sigma=1
%
%  threshold - if is different from 0, edges that are smaller than the
%  threshold are set to zero.
% 
% Output
%
%  W - is symmetrized, but no normalization is applied.
%  D - contains the row sums of W for future normalization.
%
%--------------------------------------------------------------------
% Yoel Shkolnisky 31/08/2006, Nir Sharon 2013
%--------------------------------------------------------------------

if nargin<4
    threshold=0;
end

if nargin<3
    sigma=1;
    fprintf('Using default value sigma=%d\n',sigma);
end

[NN,n] = size(distance); % in case where NN> number of points

if numel(sigma)>2
    sigma=repmat(sigma,1,NN);
end

S=exp(-distance.^2./sigma);
S=S(:);

if threshold>0
    S(S<threshold)=0;
end

I=zeros(NN,n);
for k=1:NN
    I(k,:)=1:n;
end
I=I(:);

J=index;
J=J(:);
W=sparse(I,J,S(:),n,n,n*NN);
W=(W+W.')/2;

D=sum(W,2);
end
