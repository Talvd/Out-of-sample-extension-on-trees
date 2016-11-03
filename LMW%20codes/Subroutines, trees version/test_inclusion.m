function [ arr ] = test_inclusion( T )

%test_inclusion

bfs_node_list = T.bfsiterator;
n = length(bfs_node_list);

arr = ones(n,2);
% main loop
for j=1:n
    
    % relevant ids
    parentID = bfs_node_list(j); % the current
    children = T.getchildren(parentID);
    
    if isempty(children)
        arr(j,1) = 0; 
        arr(j,2) = 0;
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
             
        % add to the its "brothers" vectors
        step = size(local_V_c,2);
        ends = ends + step;
        V_c(common_ind,starts:ends) = local_V_c;
        starts = ends+1;
        
    end
    
    arr(j,1) = OrthoCheck(V_c');
    if arr(j,1)>1e-13
        imagesc(V_c);
    end
    
    arr(j,2) = norm(V_p-V_c*V_c'*V_p);   % is V_p \subset  V_c ?
    
    if arr(j,2)>1e-13
        imagesc(V_c);
    end
 
end

end

