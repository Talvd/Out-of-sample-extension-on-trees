function [ index, distance ] = BruteForceNN( X,NN )
%
% Call the c implementation for brute force
%
[index,distance]=bfsnn2(X,NN);


index = index + 1;
distance = sqrt(distance);

end

