function [ ind, dis ] = Correlation_NN( X, NN )
% we arap the Matlab k-NN procedure with correlation metric
%
% Nir Sharon, June 2013

[ ind, dis ]  = knnsearch(X',X','K',NN,'Distance','correlation');

% an alternative - using costum made function sush as "Imagemetric"
%[ ind, dis ]  = knnsearch(X',X','K',NN,'Distance',@Imagemetric);


end

