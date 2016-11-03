function [ U_basis ] = ABCR_Version4_Overall( S,k,l )
% Main ABCR implementation - Full version
% S = data array of length n = k*2^l
% First step is done seperatly
% We have to keep the U,M values for the reccursion computition,namely
% U(::,l,i)==U(:,:,level,index)
% The over all size of input,
n = k*2^l;
% --------------- Data Validation ------------------------------
if (n ~=length(S))
    x = length(S);
    disp(['Wrong data inputs :', num2str(x),'<>', num2str(k),'*2^',num2str(l) ]);
    U_basis=0;
    return;
end
% --------------- End Data Validation ---------------------------
% --------------- Initialize step -------------------------------
U_main = zeros(n,n,l);
U_l=zeros(k,2*k,l,(n/(2*k)));
U_u=zeros(k,2*k,l,(n/(2*k)));
M=zeros(2*k,2*k,l,(n/(2*k)));
UM = zeros(2*k,2*k);  % utility array
% --------------- End Initialize step ----------------------------
% --------------- First step ("1"+"2") ---------------------------
for i=1:(n/(2*k))
    s_i=2*k*(i-1);
    muu = mu(S,k,1,i);       %(S(1+s_i)+S(i*k*2))/2;
    sigmaa = sigma(S,k,1,i); %(S(i*k*2)-S(1+s_i))/2;
    % Construct M, first level
    for deg=0:(2*k-1)  
        UM(:,deg+1) = ((S(s_i+1:s_i+2*k)-muu)/sigmaa).^(deg);
    end
    M(:,:,1,i)=UM(:,:);
    U = Ortho(UM);
    U = U';
    [U_u(:,:,1,i),U_l(:,:,1,i)] = V_partition(U);
end
% Construct the first level array
Upper = U_l(:,:,1,1);
Lower = U_u(:,:,1,1);
% Concatenation to block diagonal matrix
for i=2:(n/(2*k))
    Upper = blkdiag(Upper,U_l(:,:,1,i));
    Lower = blkdiag(Lower,U_u(:,:,1,1));
end
U_main(:,:,1) = [Upper ; Lower];
% --------------- Second step ("3") ------------------------------
% INITIALIZE.
upper_part = zeros(k,2*k);
lower_part = zeros(k,2*k);
S_1 = zeros(2*k,2*k);
S_2 = zeros(2*k,2*k);
UM = zeros(2*k,2*k);  % utility array
% INITIALIZE ENDS.
% Main loop:
for j=2:l
    for i=1:(n/((2^j)*k))
        upper_part = U_u(:,:,j-1,2*i-1)*M(:,:,j-1,2*i-1);
        lower_part = U_u(:,:,j-1,2*i)*M(:,:,j-1,2*i);
        % The S-matrices construction
        a = (mu(S,k,j,i)-mu(S,k,j-1,2*i-1))/sigma(S,k,j-1,2*i-1) ;
        b = sigma(S,k,j,i)/sigma(S,k,j-1,2*i-1);
        S_1 = S_matrix(k,a,b);
        a = (mu(S,k,j,i)-mu(S,k,j-1,2*i))/sigma(S,k,j-1,2*i) ;
        b = sigma(S,k,j,i)/sigma(S,k,j-1,2*i);
        S_2 = S_matrix(k,a,b);
        % Multiply to have M'j,i
        UM = [ upper_part*S_1 ; lower_part*S_2 ];
        M(:,:,j,i)= UM;
        U = Ortho(UM);
        U = U';
        [U_u(:,:,j,i),U_l(:,:,j,i)] = V_partition(U);
    end
    % Construct the j-level array
    Upper = U_l(:,:,j,1);
    Lower = U_u(:,:,j,1);
    % Concatenation to block diagonal matrix
    for i=2:(n/((2^j)*k))
        Upper = blkdiag(Upper,U_l(:,:,j,i));
        Lower = blkdiag(Lower,U_u(:,:,j,1));
    end
    left = n*(1-1/(2^(j-1)));
    U_main(:,:,j) = blkdiag(eye(left),[Upper ; Lower]);
end
%------- Final Basis Matrix -----------------------------
U_basis = eye(n);
for j=1:l
    U_basis=U_main(:,:,j)*U_basis;
end

end

