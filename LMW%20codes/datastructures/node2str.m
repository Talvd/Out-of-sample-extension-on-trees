% An example function to convert a tree key to a string.
% To be used with the function plot1 in class tree.
%
% Yoel Shkolnisky, June 2013.

function str=node2str(key)
I=key{:}; l=I(1); h=I(2);
str=sprintf('[%4.2f,%4.2f]',l,h);
