% Ternary tree - demo for the classes tree and queue.
% 
% Yoel Shkolnisky, June 2013.

Q=queue(10000);
depth=2;
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
    
    
    dx=(h-l)/3;
    leftinterval=[l,l+dx];
    midinterval=[l+dx l+2*dx];
    rightinterval=[l+2*dx h];
    
    [T,lnode]=addnode(T,nodeID,leftinterval);
    [T,mnode]=addnode(T,nodeID,midinterval);
    [T,rnode]=addnode(T,nodeID,rightinterval);
    
    Q=add(Q,lnode);
    Q=add(Q,mnode);
    Q=add(Q,rnode);
end
    
plot1(T,@node2str);