function [W,D]=Graph_Matrix_Var(index,distance,d)
%
% Create a symmetric sparse graph matrix where
%
%   W(i,j)=exp(-norm(X(i)-X(j))^2/sigma)
%
% Input
%
%  index,distance - are the info for the nearest neighbors of each point.
% 
% Output
%
%  W - is symmetrized, but no normalization is applied.
%  D - contains the row sums of W for future normalization.
%
%--------------------------------------------------------------------
% Yoel Shkolnisky 31/08/2006, Nir Sharon 2013, Tal Van Dijk 2015
%--------------------------------------------------------------------

[NN,n] = size(distance); % in case where NN> number of points

p0 = (sum(distance.^2)/(NN-1)).^(.5);

q0 = zeros(NN,n);
for i=1:NN
   q0(i,:)=exp(-distance(i,:).^2./2.*p0.*p0(index(i,:))); 
end

q = (2*pi)^(-d*.5)*sum(q0)./((p0.^d)*NN);
p = q.^(-.5);


% p = (sum(distance.^2)/(NN-1)).^(.5);

S = zeros(NN,n);
for i=1:NN
   S(i,:)=exp(-distance(i,:).^2./2.*p.*p(index(i,:))); 
end
S=S(:);

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
