function U = nystrom_extension(points, eps, p, U, V, flag)

%---------------------------------------------------------------
%input:
% p - the new points
% eps - epsilon to claculate the kernel
% U - the eigenvector we want to extend
% V - appropriate eigenvalue
% flag - if 1, we need to normalized the kernal. 0 - not normalization
%---------------------------------------------------------------

    if nargin<6
        flag = 1;
    end

   [m,n] = size(U);

    for l=1:m
        Ker(l) =  exp(-(norm(points(:,l)-p,2)^2)/eps);
    end
    if (flag ==1)
        Ker = Ker/sum(Ker);
    end
    for j=1:n
        Nsum =0;      
        for k=1:m
            Nsum = Nsum +  Ker(k)*U(k,j);
        end;
        U(m+1,j) = (1/V(j,j))*Nsum;
    end
end