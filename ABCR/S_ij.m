function [ Sij ] = S_ij(i,j, mu,sigma )
% The shifting and scaling operator
if (i>j)
    Sij=0;
else
    c = nchoosek(j-1,i-1);
    Sij = c*(  ((-1*mu)^(j-i))/((sigma)^(j-1)) );
end

end

