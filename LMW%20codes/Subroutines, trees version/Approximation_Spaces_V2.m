function [ T ] = Approximation_Spaces_V2( points, k, T, NN_func, NN, func_eps, metric)
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
    func_eps=[];
end

if nargin<7
    metric=[];
end

% get the list of nodes
bfs_node_list = T.bfsiterator;
n = length(bfs_node_list);

% the first app. space, that is V_0
rootID = bfs_node_list(1);
index = T.getkey(rootID);
V = GLE(points(:,index), NN_func, NN, k, func_eps, metric);
[V,~,~] = svd(V,0);
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
    c_ind  = T.getkey(currentID);  %in this point, the cell should only have the indices set
    
    % how many vectors are needed
    if isempty(T.getchildren(currentID))
        needed = length(c_ind); % if the node a leaf
    else
        needed = k;             % otherwise
    end
    
    % the restriction
    [~, rest_ind] = intersect(p_ind,c_ind);
    p_restriction = p_cell{2}(rest_ind,:);
    
    % raising the condition number
    [u_rest, s, ~]  = svd(p_restriction,0);
    u_rest          = u_rest(:,(diag(s)>eps));
    restricted_rank = size(u_rest,2);
    
    if mod(j,200)==0
        disp(num2str(j));
    end

    if j>5264
        flaggos =1;
    end
    
    % additional local vectors, if needed
    if restricted_rank<needed
        
        % local "information"
        local = GLE(points(:,c_ind), NN_func, NN, needed, func_eps, metric);
        
%         % we must pick the restricted vectors first
%         if size(p_restricted_V,2)~=restricted_rank
%             [~, ~, e] = qr(p_restricted_V,0);
%             p_restricted_V = p_restricted_V(:,e(1:restricted_rank));
% 
%             % ?? needed (01-07) ??
%             normalize = (sum(p_restricted_V.*p_restricted_V)).^(-.5);
%             p_restricted_V = p_restricted_V*diag(normalize);
%             % ?? needed (01-07) ??
% 
%         end
        
        place = needed-restricted_rank +1;
        merge = [u_rest, local(:,2:place) ];
        
        [merge, s, ~] = svd(merge,0);
        merge         = merge(:,(diag(s)>eps));
        merge_rank    = size(merge,2);
        
        while merge_rank<needed
            ind_add = needed-merge_rank;
            add     = ((place+1):(ind_add + place));
            merge   = [merge, local(:,add)];
            place   = place + ind_add;
            
            [merge, s, ~] = svd(merge,0);
            merge         = merge(:,(diag(s)>eps));
            merge_rank    = size(merge,2);
        end
        
%         %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%
%         if size(merge,2)~=merge_rank  % keep a basis
%             [~, ~, e] = qr(merge,0);
%             merge = merge(:,e(1:merge_rank));
% 
%             % ?? needed (01-07) ??
%             normalize = (sum(merge.*merge)).^(-.5);
%             merge = merge*diag(normalize);
%         end
%         %%%%%%%%%%% 01-07 %%%%%%%%%%%%%%%

        c_val = merge;

    else
        c_val = u_rest;
    end
    
    
%     % DEBUGIING ================
%     [q, ~, ~] = svd(c_val,0);
%     
%     if norm(P_rest-q*q'*P_rest)>20*eps
%         str = ['The node ',int2str(j),' with ', num2str(norm(P_rest-q*q'*P_rest)) ];
%         disp(str);
%     end
%         
%     % DEBUGIING ================
    
    
    % update the tree
    val = {c_ind, c_val};
    T = setkey(T,currentID,val);
end

end

