function [ W ] = local_ortho_V2( V , detail_V )
%---------------------------------------------------------------
% calculate W s.t.
%           W + V = detail_V,  W'*V=0

% based on iterated (classic) Gram Schmidt - orthogonelize the columns of Q to A
% A = QR, QQ' = I
%
%---------------------------------------------------------------
% Nir Sharon, June 11.
%---------------------------------------------------------------

[n,m] = size(detail_V);

if size(V,1)~=n
    error('Wrong dimensions in direct differences operation')
end

% initializing
IterLimit = 10;
needed    = size(detail_V,2)-size(V,2);
w_ind     = 1;
W         = zeros(n,needed);
e_ind     = 1;
extra     = zeros(n,needed);
Q1        = V;

% main loop
for j=1:m
    % the current vector to handle
    q0 = detail_V(:,j);
    
    % first iteration
    s0 = Q1'*q0;
    q = q0 - Q1*s0;     % project q0 into V^\bot
    
    % ensure q in the space (17-06)
    q = detail_V*(detail_V'*q);  
    % --------------------- (17-06)
        
    % refine the orthogonality (iterated GS)
    count = 1;
    while ((norm(Q1'*q)>eps)&&(count<IterLimit))
        s0 = Q1'*q;
        q = q - Q1*s0;
        count = count+1;
    end
    
    % New values and normalization
    if (norm(q)>eps)&&(norm(Q1'*q)<=eps)
        W(:,w_ind) = q * (1/(norm(q)));
        w_ind = w_ind + 1;
        % The relevant vectors which are not perfectly orthogonal
    else if (norm(q)>eps)
            extra(:,e_ind) = q * (1/(norm(q)));
        end
    end
end %main loop

% final setting of W
w_ind = w_ind-1;       % number of actual vectors in W
W     = W(:,1:w_ind);

% if we picked "too many" vectors
if w_ind>needed
    [W,dd,~] = svd(W,0);
    W = W(:,diag(dd)>eps);
    if size(W,2)>needed
        W = W(:,1:needed);
    end
end

% if we picked "too few" vectors
if size(W,2)<needed
    
    w_ind = size(W,2);
    add = needed-w_ind;
    
    % adding from "extra" array
    if size(extra,2)<add
        add = size(extra,2);
        W(:,(w_ind+1):(w_ind+add)) = extra;
        warning('Failure of local orthogonality-not enough orthogonal vectors');
    else  
        ortho_msr = V'*extra;
        [~, r_ind] = sort(ortho_msr);
        % adding the most orthogonal vectors available
        W(:,(w_ind+1):(w_ind+add)) = extra(:,r_ind(1:add));
    end
end


end

% % %test section
% %
% % if (size(W,2)+size(V,2))~=size(detail_V,2)
% %     warning('severe rank issues');
% % end
% %
% % orthog = norm(V'*W);
% % if orthog>(100*eps)
% %     warning('ortho accuarcy of %d', orthog);
% % end


