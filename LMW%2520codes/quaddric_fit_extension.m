function [V] = quaddric_fit_extension(points,extra_p,V,NN, flag)

if nargin<4
    NN = 10;
end

if nargin<5
    flag = 0;
end

d = size(points,1);
  
[n, m] = size(V);
if (m ==1 && flag ~= 1)
    V(n+1,1)=V(n,1);
else    
    if (NN > n)
        NN = n;
    end
    idx= knnsearch(points',extra_p','k',NN);

    for i=1:size(V,2)       
        vi = V(1:n,i);
        if (d==1)
            [ A, a, b ] = quaddric_fit_lscov(points(:,idx), vi(idx) ,extra_p);
        else
            [ A, a, b ] = quaddric_fit( points(:,idx), vi(idx));
        end
       
        V(n+1,i) = extra_p'*A*extra_p + a'*extra_p + b;
    end
end

end

