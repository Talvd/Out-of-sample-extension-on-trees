function [ A, a, b ] = quaddric_fit(x, f)
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

number_monom = nchoosek(r_d+2,2) - r_d;

C = zeros(n,number_monom);
A = zeros(r_d);
% Constructing the system
for k=1:n
      C(k,1:1+r_d)  = [1 x(:,k)'];
    % quadric coef
    count=2+r_d;
    for i=1:r_d
        for j=i:r_d          
            if j>i
                C(k,count) =  x(i,k)*x(j,k);
                count = count+1;
            end            
        end
    end 
end

q = regress(f,C);

a = q(2:1+r_d);

b = q(1);

count=2+r_d;
for i=1:r_d
    for j=i:r_d       
        if j>i
            A(i,j) = .5*q(count);
            A(j,i) = .5*q(count);
            count = count+1;
        end
       
    end
end    
end

