function [ T ] = Approximation_Spaces( points, k, T, NN_func, NN, metric)
%---------------------------------------------------------------
% The first step of the algorithm - constructing the nested spaces
%
% Input
%  points - the data points: d X N where d is the dimension
%  k      - the number of minimal points per partition subset.
%  T      - the heirarchy tree
%  NN     - Number of nearest neighbours to consider in the affinity
%           construction matrix
%
% Output
%  T - a tree class variable
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

if nargin<6
    metric=[];
end

% get the list of nodes
bfs_node_list = T.bfsiterator;
n = length(bfs_node_list);

% the first app. space, that is V_0
rootID = bfs_node_list(1);
index = T.getkey(rootID);
V = GLE(points(:,index), NN_func, NN, k, metric);
val = {index, V};
T = setkey(T,rootID,val);

% rest of nodes follow
for j=2:n
    
    % relevant ids
    currentID = bfs_node_list(j);
    parentID  = T.getparent(currentID);
    
    % retrieve information
    p_cell = T.getkey(parentID);
    p_ind  = p_cell{1};
    c_ind  = T.getkey(currentID);
    
    % how many vectors are needed
    if isempty(T.getchildren(currentID))
        needed = length(c_ind); % if the node a leaf
    else
        needed = k;             % otherwise
    end
    
    % the restriction
    [~, rest_ind] = intersect(p_ind,c_ind);
    p_restricted_V = p_cell{2}(rest_ind,:);
    
    P_rest = p_restricted_V;
    
    %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%
    p_restricted_V = p_restricted_V*diag((sum(p_restricted_V.*p_restricted_V)).^(-.5));  %30-06       
    %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%

    % additional local vectors, if needed
    restricted_rank = rank(p_restricted_V,eps);
    if restricted_rank<needed
        local = GLE(points(:,c_ind), NN_func, NN, needed, metric);
        
        % we must pick the restricted vectors first
        if size(p_restricted_V,2)~=restricted_rank
            [~, ~, e] = qr(p_restricted_V,0);
            p_restricted_V = p_restricted_V(:,e(1:restricted_rank));

            % ?? needed (01-07) ??
            normalize = (sum(p_restricted_V.*p_restricted_V)).^(-.5);
            p_restricted_V = p_restricted_V*diag(normalize);
            % ?? needed (01-07) ??

        end
        
        place = needed-restricted_rank +1;
        merge = [p_restricted_V, local(:,2:place) ];
        merge_rank = rank(merge,eps);
        while merge_rank<needed
            add = ((place+1):(needed-merge_rank + place));
            merge = [merge, local(:,add)];
            merge_rank = rank(merge,eps);
        end
        
        %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%
        if size(merge,2)~=merge_rank  % keep a basis
            [~, ~, e] = qr(merge,0);
            merge = merge(:,e(1:merge_rank));

            % ?? needed (01-07) ??
            normalize = (sum(merge.*merge)).^(-.5);
            merge = merge*diag(normalize);
        end
        %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%

        c_val = merge;

    else
        c_val = p_restricted_V;
    end
    
    
    % DEBUGIING ================
    [q, ~, ~] = svd(c_val,0);
    
    if norm(P_rest-q*q'*P_rest)>20*eps
        str = ['The node ',int2str(j),' with ', num2str(norm(P_rest-q*q'*P_rest)) ];
        disp(str);
    end
        
    % DEBUGIING ================
    
    
    % update the tree
    val = {c_ind, c_val};
    T = setkey(T,currentID,val);
end

end

