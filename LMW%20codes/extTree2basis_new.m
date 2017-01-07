function [ S ] = extTree2basis_new( T, branch, n)
%---------------------------------------------------------------
% The final step of the algorithm - extract the elements of the new basis row from
% the tree.
% 
% Input
%  T - a tree class variable. The data (key) of each node is a 1X2 or 1X3 
%      cell, according to being a leaf or not, respectively.
% branch - the index of the relevent nodes
% n - the length of the row
%
% Output
%  S - the new basis row, pack in sparse matrix
%
%---------------------------------------------------------------
% Tal Van Dijk Dec 2016
%--------------------------------------------------------------

C = zeros(1,n);
bs = size(branch,2);

% the only space from the nested spaces 
val = T.data{1};
[~, k] = size(val{2});
C(:,1:k) = val{2}(end,:);

for i=bs:-1:2 
    val = T.data{branch(i)};
    [~,ik] = size(val{3});
    try
        curr_id= val{4};
        C(:, curr_id : curr_id + ik -1) = val{3}(end,:); 
    catch ex
        jj=1;
    end
end

S = sparse(C);

end

