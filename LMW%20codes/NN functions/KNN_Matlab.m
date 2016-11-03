function [ ind, dist ] = KNN_Matlab( data, NN)
% wrapping the matlabs function
% data -  d X N

[ind, dist ] = knnsearch(data',data','K',NN);

ind  = ind';
dist = dist';

end

