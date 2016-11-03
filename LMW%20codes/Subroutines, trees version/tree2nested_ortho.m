function [ C ] = tree2nested_ortho( T )
%---------------------------------------------------------------

%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------


% get the list of nodes
bfs_node_list = T.bfsiterator;
m = length(bfs_node_list);


val = T.getkey(bfs_node_list(1));

[N,k] = size(val{2});

% main loop
C = [];
for j=1:m
    
    val = T.getkey(bfs_node_list(j));
    ind = val{1};
    V = zeros(N,k);
    [V(ind,1:(size(val{2},2))), ~,~] = svd(val{2},0); 
        C = [C,V];
end


end

