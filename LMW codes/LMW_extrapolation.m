function [ B, T ] = LMW_extrapolation( points,p, T,B, NN, extension_func)
%---------------------------------------------------------------
% Extend LMW basis for the new point p.
%
% Input
%  points - the data points: d X N where d is the dimension
%  p - the new point.
%  k - the number of minimal points per partition subset.
%  T - the tree LMW basis  
%
% Output
%  B - the LMW basis
%

if (nargin<6)||isempty(extension_func)
    extension_func = @quaddric_fit_extension;
end

n = size(points,2);

branch = find_branch( p,points, T, NN);

%extend root
bs = size(branch,2);
parentID=branch(bs);    
p_cell = T.getkey(parentID);
idx = 1:n;
k = size(p_cell{2},2);
if (k > 1)
    idx= knnsearch(points',p','k',NN);
end
V = p_cell{2};
V = extension_func(points(:,idx),p,V(idx,:));
p_cell{2} = [p_cell{2};V];

p_cell{1} = [p_cell{1} n+1];
T = setkey(T,parentID,p_cell); 
for i=bs-1:-1:1 
    currentID=branch(i);
    c_cell = T.getkey(currentID);
    
    V =  c_cell{2};
%     NN = 10;
%     if (i==1)
%         NN = 3;
%     end

    curr_p=points(:,c_cell{1});
    if ( k>1)
        [~,curr_idx] = intersect(c_cell{1}, idx);
        if (size(curr_idx,2) > size(curr_p,2))
            curr_idx = 1:size(curr_p,2);
        end
    else
        curr_idx = 1:size(curr_p,2);
    end
    if(isempty(curr_idx))
        curr_idx= knnsearch(curr_p',p','k',NN);
    end
    V= extension_func(curr_p(:,curr_idx),p,V(curr_idx,:));
    c_cell{2} = [c_cell{2}; V];
    c_cell{1} = [c_cell{1} n+1];
    T = setkey(T,currentID,c_cell); 
  
    if(size(p_cell,2) >= 3)       
        child_ind = T.getchildren(parentID);
        V = [];
        for j=child_ind
            child = T.getkey(j);
            Vj = zeros(n+1,size(child{2},2));
            Vj(child{1},:)=child{2};
            V = [V Vj];
        end
        w = zeros(n, size(p_cell{3},2));        
        w(p_cell{1}(1:end-1),:) = p_cell{3};
        v =V(1:end-1,:);       
%         [v ,z, u] = svd(v,0);
% 
%         r= v'* w;  % coeffiecients of the orthogonal representation
%         V = V*u*inv(z);
%         W = V*r; 
        b= v'*w;
        W=V*b;
        p_cell{3} = W(p_cell{1},:);
        T = setkey(T,parentID,p_cell); 
        
        p_cell = c_cell;
        parentID = currentID;
    end    
end
 
% packing the basis into a sparse matrix
B_ext = extTree2basis_new(T,branch,n);
% B_ext = extTree2basis(T);

B = sparse([B;B_ext]);

end

