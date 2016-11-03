function [ V_u,V_l ] = V_partition( V )
% divide the upper and lower parts of V
n = size(V,1);
if mod(n,2)==0
    V_u = V(1:(n/2),:);
    V_l = V(((n/2)+1):n,:);
else
    V_u = V(1:((n-1)/2),:);
    V_l = V((((n+1)/2)):n,:);
end    

end

