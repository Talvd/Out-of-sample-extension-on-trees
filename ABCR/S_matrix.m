function [ S ] = S_matrix(k, mu, sigma )
% Scale and shifting matrix of 2k X 2k
S = zeros(2*k,2*k);
for i=1:2*k
    for j=i:2*k
        S(i,j)=S_ij(i,j, mu,sigma);
    end
end

end

