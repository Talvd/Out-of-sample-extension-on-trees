function [ A, a, b ] = quaddric_fit_lscov( x, f ,p)
% Given data points on parametric hyperplane
% dat, and the values of the manifold, we look for 
% best fit of the form 
%   x'Ax + a'x + b
% Input
% x - n points on the refrence (r_d X n)
% f - n values of the manifold (n X 1)
%
%-------------------------------------------------------------
% Nir Sharon, March 2013
%-------------------------------------------------------------

[r_d,n] = size(x);

number_monom = nchoosek(r_d+2,2);% r_d*(r_d+1)/2;

C = zeros(n,number_monom);
A = zeros(r_d);

% Constructing the system
for k=1:n
    % quadric coef
    count=1;
    for i=1:r_d
        for j=i:r_d
            C(k,count) = x(i,k)*x(j,k);
            if j~=i
                C(k,count) = 2*C(k,count);
            end 
            count = count+1;
        end
    end
    % the linear part
    for i=1:r_d
        C(k,end-1-r_d+i) = x(i,k);
    end
    % constant coef
    C(k,end) = 1; 
end

eps = 0.01;

for l=1:n
    d(l) =  exp(-(norm(x(:,l)-p,2)^2)/eps^2);
end
if (size(d,2) > 1)
    q = lscov(C,f,d);
else
     q = pinv(C)*f;
end


count=1;
for i=1:r_d
    for j=i:r_d
        if i==j
            A(i,j) = q(count);
        else
            A(i,j) = .5*q(count);
            A(j,i) = .5*q(count);
        end
        count = count+1;
    end
end
a = q((number_monom-r_d):(number_monom-1));

b = q(end);
   
end

