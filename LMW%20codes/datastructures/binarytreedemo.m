% Binary tree - demo for the classes tree and queue.
% 
% Yoel Shkolnisky, June 2013.

Q=queue(20);
depth=3;
T=tree;

[T,nodeID]=addnode(T,0,[0,1]);
Q=Q.add(nodeID);

while ~isempty(Q)
    [Q,nodeID]=remove(Q);
    
    key=getkey(T,nodeID);
    I=key{:};
    l=I(1);
    h=I(2);
    dd=getdepth(T,nodeID);
    
    fprintf('Procedding interval [%4.2f,%4.2f] \t depth=%d\n',l,h,dd);
        
    if dd>depth-1
        continue
    end
    
    
    leftinterval=[l,(l+h)/2];
    rightinterval=[(l+h)/2,h];
    
    [T,lnode]=addnode(T,nodeID,leftinterval);
    [T,rnode]=addnode(T,nodeID,rightinterval);
    
    Q=add(Q,lnode);
    Q=add(Q,rnode);
end
    
plot1(T,@node2str);