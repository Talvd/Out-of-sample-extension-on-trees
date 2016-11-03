% A function to convert a tree key to a string.
% To be used with the function plot1 in class tree.
%
% Nir Sharon, June 2013.

function str=node2digitmean(key)
%val = length(key); 


%labels =  [ ones(200,1) ; 8*ones(200,1)];
struct = load('partial_labels');

labels = struct.labels;


str = sprintf('#1=%d\t #8=%d',sum(labels(key)==1), sum(labels(key)==8));
%str=sprintf('%4.2f',val);
end
