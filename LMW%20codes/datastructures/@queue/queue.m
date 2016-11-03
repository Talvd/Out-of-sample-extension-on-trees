classdef queue
    % Queue data structure. Stores only numbers. Implemeted as a cyclic
    % array.
    %
    % Todo: function peek to query the head of the queue without removing
    % the head element.
    %
    % Yoel Shkolnisky, June 2013.
    
    properties (SetAccess = 'private')
        maxitems % Maximal number of nodes to be stored. Whenever this 
                 % number is reached, the size of the tree doubles. 
                        
        data % The key (integer) of each node.                
        
        head % Head of the queue
        
        tail % Tail of the queue
        
        count % Number of items in the queue
        
    end
    
    methods
        function Q=queue(N)
            % Create a queue that can hold N items.

            if N==0
                error('Cannot create a queue of size zero');
            end
            Q.data=zeros(N,1);
            Q.head=1;
            Q.tail=0;
            Q.count=0;
            Q.maxitems=N;
        end
        
        function Q2=realloc(Q,newsize)
            if Q.count>=newsize
                error('newsize is too small');
            end
            Q2=queue(newsize);
            while ~isempty(Q)
                [Q,key]=remove(Q);
                Q2=add(Q2,key);
            end
        end
            
        
        function Q=add(Q,key)
            % Add key to the queue Q.
            if Q.count>=Q.maxitems
                newsize=2*Q.maxitems;
                fprintf('Queue is full. Reallocating from size %d to size %d.\n',Q.maxitems,newsize);
                Q=realloc(Q,newsize);
            end
            Q.tail=Q.tail+1;
            Q.tail=mod((Q.tail-1),Q.maxitems)+1;
            Q.data(Q.tail)=key;
            Q.count=Q.count+1;
            %fprintf('head=%d\ttail=%d\n',Q.head,Q.tail);
        end
        
        function [Q,key]=remove(Q)
            % Remove the head of the queue Q.
            if Q.count==0
                error('Queue is empty');
            end
            key=Q.data(Q.head);
            Q.head=Q.head+1;
            Q.head=mod((Q.head-1),Q.maxitems)+1;
            Q.count=Q.count-1;
            %fprintf('head=%d\ttail=%d\n',Q.head,Q.tail);
        end
        
        function TF=isempty(Q)
            % Is the queue empty?
            
            TF=(Q.count==0);            
        end
        
        function tostring(Q)
            % Print the queue.
            h=Q.head;
            for k=1:Q.count
                fprintf('%d\t',Q.data(h));
                h=h+1;
                h=mod((h-1),Q.maxitems)+1;
            end
            fprintf('\n');
        end
    end
end