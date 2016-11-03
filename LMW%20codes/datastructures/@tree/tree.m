classdef tree
    % Tree data structure.
    % Currently supports only addition and query. key can be anything.
    % 
    % Implementation:
    % Each node in the tree has a row in the (sparse) children array which
    % contains the children of the node. The ChildrenCount contain the
    % number of children. Each node also has a pointer to its parent in the
    % Parents array.
    %
    % Yoel Shkolnisky, June 2013.
    
    properties (SetAccess = 'private')
        maxnodes % Maximal number of nodes to be stored. Whenever this 
                 % number is reached, the size of the tree doubles. 
        
        maxchildren % Maximal number of children per node.
        
        Parents  % The parent of each node.
        
        Children % A list of chidlren for each node.
        
        ChildrenCount % Number of children for each node.
        
        data % The key of each node.
        
        nnodes % Number of nodes in the tree.
        
    end
    
    methods
        
        % CONSTRUCTOR        
        function root=tree(maxnodes,maxchildren)
            % maxnodes  Maximal number of nodes in the tree Default 1000.
            % maxchildren   Maximal number of children per node. Default 3.            
            root.maxnodes=1000;
            if nargin>0
                root.maxnodes=maxnodes;
            end
            
            root.maxchildren=3;
            if nargin>1
                root.maxchildren=maxchildren;
            end
            root.Parents=zeros(root.maxnodes,1);
            root.Children=sparse(root.maxnodes,root.maxchildren);
            root.ChildrenCount=zeros(root.maxnodes,1);
            root.data=cell(root.maxnodes,1); 
            root.nnodes=0;
        end
        
        function T2=realloc(T,newsize)
            % Realloc the tree to a new size.
            
            % Check if newsize if sufficient.
            if newsize<=T.nnodes
                error('newsize is too small');
            end
            T2=tree(newsize,T.maxchildren); % Allocate the new tree.
            
            % Now copy the data.
            T2.Parents(1:T.nnodes)=T.Parents(1:T.nnodes);
            
            idx=find(T.Children);
            [I,J]=ind2sub([T.maxnodes T.maxchildren],idx);
            idx2=sub2ind([newsize T.maxchildren],I,J);
            T2.Children(idx2)=T.Children(idx);            
            T2.ChildrenCount(1:T.nnodes)=T.ChildrenCount(1:T.nnodes);
            T2.data(1:T.nnodes)=T.data(1:T.nnodes);
            T2.nnodes=T.nnodes;
        end
            
            
        
        function [T nodeID]=addnode(T,parent,key)
            % Use parent = 0 to add the root.
            if T.nnodes==T.maxnodes
                newsize=2*T.maxnodes;
                fprintf('Reached size limit %d. Reallocate tree to size %d.\n',T.maxnodes,newsize)
                T=realloc(T,newsize);
            end
            
            if parent>T.nnodes
                error('Parent does not exist');
            end
            
            T.nnodes=T.nnodes+1;
            nodeID=T.nnodes;
            
            if parent==0  % Adding the root. No parent to update.
                T.data{nodeID}=key;
            else    
                T.Parents(nodeID)=parent;
                T.ChildrenCount(parent)=T.ChildrenCount(parent)+1;
                T.Children(parent,T.ChildrenCount(parent))= nodeID;
                T.data{nodeID}=key;
            end                
            
        end
        
        function T=setkey(T,nodeID,key)
            % Set the key associated with nodeID.
            T.data{nodeID}=key;
        end
        
        function key=getkey(T,nodeID)
            % Get the key assoicates with nodeID.
            key=T.data{nodeID};
        end
        
        function parent=getparent(T,nodeID)
            % Get the parent id of the given node.
            parent=T.Parents(nodeID);
        end
        
        function children=getchildren(T,nodeID)
            % Get the list of children of nodeID.
            nc=nnz(T.Children(nodeID,:));
            children=full(T.Children(nodeID,:));
            children=children(1:nc);
        end
        
        function depth=getdepth(T,nodeID)
            % Get the depth of nodeID.
            depth=0;
            parent=getparent(T,nodeID);
            while parent~=0
                depth=depth+1;
                parent=getparent(T,parent);
            end
        end
   
        
        function nodelist=bfsiterator(T)
            % Return the list of nodes in the tree in BFS order. The
            % returned nodelist contains only the indices of the nodes.
            
            nodelist=[];
            if T.nnodes==0
                return;
            end
            Q=queue(T.nnodes);
            nodelist=zeros(T.nnodes,1);
            Q=add(Q,1);
            idx=1;
            while ~isempty(Q)
                [Q,nodeID]=remove(Q);
                nodelist(idx)=nodeID;
                idx=idx+1;
                ch=getchildren(T,nodeID);
                for k=1:numel(ch)
                    Q=add(Q,ch(k));
                end
            end
        end
        
        function plot1(T,node2str)
            % Plot the tree. 
            % The parameter node2str is a function pointer with the
            % signature str=node2str(key) that converts a key to a string.
            
            % Find how many nodes are in each depth
            nodelist=bfsiterator(T);
            maxdepth=getdepth(T,nodelist(end));
            depthcount=zeros(maxdepth+1,1);
            for k=1:numel(nodelist)
                dd=getdepth(T,nodelist(k));
                depthcount(dd+1)=depthcount(dd+1)+1;
            end
            
            % Global height and width of the figure.
            height=100*(maxdepth+1);
            width=100*max(depthcount);
            verticalinterval=100;
            
            % Populate drawing table - coordinates where to plot each node
            % according to its depth and horizontal location in its level.            
            coordmap=zeros(numel(nodelist),2); %Row nodeID contain x and y for that node. 
            currdepth=1;
            horizontalinterval=width/2;
            idx=1; % This is the drawing position of the current node in 
                   % the current depth. Indexes the nodes in the current 
                   % depth to determine where to plot each node.
            for k=1:numel(nodelist)
                if getdepth(T,nodelist(k))~=currdepth-1
                    % If we have advanced to the next level, update the
                    % horizonal spacing between nodes.
                    currdepth=currdepth+1;
                    horizontalinterval=width/(depthcount(currdepth)+1);
                    idx=1; % Depth has increased, so reset the drawing 
                           % position to the begining of the line.
                end
                x=idx*horizontalinterval;
                y=currdepth*verticalinterval;
                coordmap(nodelist(k),1)=x;
                coordmap(nodelist(k),2)=y;
                idx=idx+1;
            end
            
            %Draw nodes
            clf;
            hold on;
            axis ij
            for k=1:numel(nodelist)
                x=coordmap(k,1);
                y=coordmap(k,2);
                %scatter(x,y,50,'r','filled')
                key=getkey(T,nodelist(k));
                str=node2str(key);
                text(x,y,str,'HorizontalAlignment','center');
            end
            hold off;
            axis([0 width verticalinterval-10 height]);
            axis off
        
            % Draw lines
            hold on;
            axis ij
            for k=1:numel(nodelist)
                % Get coordinates of current node.
                x=coordmap(k,1);
                y=coordmap(k,2);
                parent=getparent(T,nodelist(k));
                if parent==0
                    continue;
                end
                
                % Get the coordinates of the parent
                px=coordmap(parent,1);
                py=coordmap(parent,2);
                
                % Vector from parent to current node.
                dx=x-px;
                dy=y-py;
                
                % Shorten the vector by 10% at each end so the lines don't
                % overlap with the text of the node.
                x0=px+0.1*dx; x1=px+0.9*dx;                
                y0=py+0.1*dy; y1=py+0.9*dy;
                
                % Plot line.
                line([x0 x1],[y0 y1]);
            end
            hold off;
                        
%             % Draw
%             clf;
%             hold on;
%             axis ij
%             verticalinterval=100;
%             for k=1:maxdepth+1
%                 horizontalinterval=300/(depthcount(k)+1);
%                 for j=1:depthcount(k)
%                     x=j*horizontalinterval;
%                     y=k*(verticalinterval-1);
%                     scatter(x,y,50,'r','filled')
%                 end
%             end
%             hold off;
%             %axis([0 300 0 (maxdepth+1)*verticalinterval]);
%             axis off
        end
    end
end