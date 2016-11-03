function [ T ] = Fast_Orthogonalization(T)
%---------------------------------------------------------------
% The second step of the algorithm - orthogonalization
% 
% Input
%  T - a tree class variable. Assume the data (key) of each node is a 1X2 cell
%
% Output
%  T - a tree class variable after the orthogonalization
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

% get the list of nodes
bfs_node_list = T.bfsiterator;
n = length(bfs_node_list);

% % % prepare the root 
% % rootID    = bfs_node_list(1);
% % root_cell = T.getkey(rootID);
% % [V_p, ~, ~]  = svd(root_cell{2},0);
% % root_cell{2} = V_p;
% % T = setkey(T,rootID,root_cell);
% % 
% % %============ DEBUGGING ONLY ============
% % if (OrthoCheck(V_p')<1e-14);  
% %     warning('Root ortho');
% % end
% % %============DEBUGGING ONLY============


% main loop
% During process: V_p = p_cell{2} is always (suppose to be) an orthogonal set.
for j=1:n
    
    % relevant ids
    parentID = bfs_node_list(j); % the current
    children = T.getchildren(parentID);
    
    if isempty(children)
        continue
    end
    
    % we construct the relevant nesting spaces
    p_cell = T.getkey(parentID);
    p_ind  = p_cell{1}; 
    V_p    = p_cell{2};
    
%     [V_p, ~, ~] = svd(V_p,0);
%     assert(OrthoCheck(V_p')<10e-14);

    
    % this matrix will include the basis of the detalied space
    V_c = zeros(numel(p_ind),1);
  
    % pass over the children to construct the "detailed space"
    starts = 1;
    ends = 0;
    for l=1:numel(children)
        
        % extract data
        c_cell    = T.getkey(children(l));
        c_ind     = c_cell{1};
        local_V_c = c_cell{2};
        
        % adapt indices
        [~, common_ind] = intersect(p_ind,c_ind);
       
        % setting the detalied space 
% %         OrthoCheck(local_V_c');
% %         [local_V_c, ~, ~] = svd(local_V_c,0);   % CHECK
        
% % DEBUGGING +++++++++++++++++++++++
%         orth = OrthoCheck(local_V_c');
%         disp(['Local orthogonality ', num2str(orth)]);
%         if orth>1e-14
%             warning('problem');
%         end
% % DEBUGGING +++++++++++++++++++++++
      
        
        % add to the its "brothers" vectors
        step = size(local_V_c,2);
        ends = ends + step;
        V_c(common_ind,starts:ends) = local_V_c;
        starts = ends+1;
        
% %         % update the child's key
% %         c_cell{2} = local_V_c;
% %         T = setkey(T,children(l),c_cell);
    end
       
%     V_p = p_cell{2};
%     [V_p, ~, ~] = svd(V_p,0);
%     assert(OrthoCheck(V_p')<10e-14);
    
    % direct substruction procedure
    W_p = local_ortho_V2(V_p,V_c);

%     % % DEBUGGING +++++++++++++++++++++++
%     disp('is W_p \bot V_p ?');
%     norm(W_p'*V_p)           % is W_p \bot V_p ?
%     disp('is W_p \in  V_c ?');
%     norm(W_p-V_c*V_c'*W_p)   % is W_p \in  V_c ?
%     % % DEBUGGING +++++++++++++++++++++++
    
    
    p_cell{3} = W_p;
    T = setkey(T,parentID,p_cell);
end

end

