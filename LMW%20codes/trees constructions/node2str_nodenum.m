% A function to convert a tree key to a string.
% To be used with the function plot1 in class tree.
%
% Nir Sharon, June 2013.

function str=node2str_nodenum(key)
val = length(key); 


% l=I(1); h=I(2);
str=sprintf('%d',val);
