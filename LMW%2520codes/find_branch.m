function [ branch ] = find_branch( p,points, T, NN)
%find branch for the new point

% Input
%  p - the data point: d X 1 where d is the dimension
%  T      - the heirarchy tree
%
% Output
%  branch - array of the node of the branch

%find the distance between the new point p and all the original points
% n = size(points,2);
% [index, d ] = knnsearch(points',p','K',n);
%     
% dist=zeros(1,n);
% dist(:,index) = d;
% Q=queue(T.nnodes);
% Q=add(Q,1);
% idx=1;
% while ~isempty(Q)
%     [Q,nodeID]=remove(Q);
%     idx=idx+1;
%     ch=T.getchildren(nodeID);
%     if(isempty(ch))
%         current_ind=T.data{nodeID}{1};
%         [~,ia,~] = intersect(ind,current_ind);
%         current_dist = dist(ia);
%         avg_dist = mean(current_dist);
%         p_dist = [p_dist; nodeID  avg_dist];
%     else        
%         for k=1:numel(ch)
%             Q=add(Q,ch(k));
%         end
%     end
% end

j=1;
ind=1;
val = T.data{j};
p_dist = Inf*ones(size(points,2),2);
while (~isempty(val))    
    l = length(val);   
    if( l==2)
%        [~,ia,~] = intersect(ind,val{1});       
%        current_dist = dist(ia);
       current_dist= sqrt(sum((points(:,val{1}) - p*ones(1,size(val{1},2))).^2,1));       
       avg_dist = mean(current_dist);
       p_dist(ind,:)= [j  avg_dist];
       ind=ind+1;
    end
    j=j+1;    
    val = T.data{j}; 
end

m= min(p_dist(:,2));
leaf_ind = p_dist(:,2)==m;
leaf = p_dist(leaf_ind,1);
leaf= leaf(1,1);

branch = leaf;
currentID = leaf;
while (currentID ~=1)
     currentID = T.getparent(currentID);
     branch =[branch currentID];     
end

end


