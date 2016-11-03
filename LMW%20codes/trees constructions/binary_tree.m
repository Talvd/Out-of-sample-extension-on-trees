function [ T ] = binary_tree(points, binary_partition, k, func_NN, NN, plotting)
%--------------------------------------------------------------- 
% Construct a binary hierarchy tree using binary partition method (binary_partition)  
% We use the tree structure
%
% Input
%  points - the data points: d X N where d is the dimension
%  k - the number of minimal points per partition subset.
%  NN - Number of nearest neighbours to consider in the affinity 
%       construction matrix
%       
% Output
%  T - a tree class variable
%
%---------------------------------------------------------------
% Nir Sharon June 2013
%---------------------------------------------------------------

if nargin<6
    plotting = 1;
end
n = size(points,2);
ind_arr = 1:n;
depth = 2*ceil(log(n/k));
Q=queue(ceil(n/k));
T=tree;

[T,nodeID]=addnode(T,0,[1:n]);
Q=Q.add(nodeID);


while ~isempty(Q)
    [Q,nodeID]=remove(Q);
    
    index=getkey(T,nodeID);
    
    ind = (index);  %updated!
    
    if iscell(ind)
        r=1;
    end

	[left_index, right_index] = binary_partition(points(:,ind),func_NN,NN);


    dd=getdepth(T,nodeID);
  
    if (numel(ind_arr(left_index))<k)||(numel(ind_arr(right_index))<k) % does only numel(left_index) enough?
        dd=depth;
    end

        %     
%     fprintf('Procedding interval [%4.2f,%4.2f] \t depth=%d\n',l,h,dd);
        
    if dd>depth-1
        continue
    end
    
    current_arr = ind_arr(ind);
    left_index = current_arr(left_index);
    right_index = current_arr(right_index);
    
    
    [T,lnode]=addnode(T,nodeID,left_index);
    [T,rnode]=addnode(T,nodeID,right_index);
    
    Q=add(Q,lnode);
    Q=add(Q,rnode);
end
if plotting
    plot1(T,@node2str_nodenum);
end

end

