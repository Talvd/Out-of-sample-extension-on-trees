function [s Acs P] = deterministicInterpolativeDecomposition(mat,k)
% deterministic interpolative decomposition algorithm
% INPUT: 
% mat - the real matrix we wish to decompose
% k - the number of columns to return. The numerical rank

% OUTPUT:
% Acs - A colums skeleton. a k-subset of the culomns of A. Also refered as
% B). 
% P - the interpulation matrix. Subset of the columns of P is the matrix
% I
% s - the chosen columns of A.

% Note: If we multiply Acs and P we get an approximation to A.  A ~ Acs*P 

[l n] = size(mat);
[q r Pr] = qr(mat); % q:l x l,    r: l x n,    Pr: n x n
Q = q(:,1:k); % l x k
r11 = r(1:k,1:k); %  k x k
r12 = r(1:k,k+1:end); % k x (n-k)
Acs = Q*r11; % l x k  - submatrix of Y
T = r11\r12; % k x (n-k)
P = [eye(k) T]*Pr.'; %P: k x n    %Y ~ ZP
s = Pr(:,1:k).'*(1:n).';