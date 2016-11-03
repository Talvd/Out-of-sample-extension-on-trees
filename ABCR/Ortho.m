function [ W ] = Ortho( V )
%Orthogonelize the columns of V to W
%m = size(V,1); % size of the columns of V
n = size(V,2); % number of columns in V
epsilon = 0.01;
%Graham Shmidt process
for j=1:n
    for t=(j-1):-1:1
        V(:,j)=V(:,j)-(V(:,j)'*V(:,t))*V(:,t);
    end
    if (norm(V(:,j))>epsilon)
        V(:,j)=(1/(norm(V(:,j))) )*V(:,j);
    end
end
W=V;

end

