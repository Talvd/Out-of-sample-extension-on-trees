function [ S ] = extTree2basis( T )
%---------------------------------------------------------------
% The final step of the algorithm - extract the elements of the basis from
% the tree.
% 
% Input
%  T - a tree class variable. The data (key) of each node is a 1X2 or 1X3 
%      cell, according to being a leaf or not, respectively.
%
% Output
%  S - the basis, pack in sparse matrix
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------


% get the list of nodes
% bfs_node_list = T.bfsiterator;
% rootID = bfs_node_list(1);
% m = length(bfs_node_list);
% 
% % the only space from the nested spaces 
% val = T.getkey(rootID);
% k = size(val{2},2);
% %[B(:,1:k), ~, ~] = svd(val{2},0); changed 17-6, already done in step 2
% B(:,1:k) = val{2};
% N = size(B,1);

% main loop
C = [];
j=1;
% the only space from the nested spaces 
val = T.data{j};
[n, k] = size(val{2});
B(1,1:k) = val{2}(n,:);

while (~isempty(val)) 
    p = length(val);   
    if( p==3)        
        ik = size(val{3},2);
        V = zeros(1,ik);
        if (val{1}(end) == n)
            V(1,1:(size(val{3},2))) = val{3}(end,:); 
        end        
        C = sparse([C,V]);
    end
    j=j+1;
    val = T.data{j}; 
end

S = sparse([B,C]);

end

