function [s B P] = randomInterpolativeDecomposition(A,k)
% random  interpolative decomposition algorithm%

% INPUT: 
% A - the real matrix we wish to decompose
% k - the number of columns to return. The numerical rank

% OUTOUT:
% B - a k-subset of the culomns of A. 
% P - the interpulation matrix. Subset of the columns of P is I
% s - the chosen columns of A.

% NOTE: If we multiply B and P we get an approximation to A.  A ~ B*P 

[m n]  = size(A);
l = min([k+3 m n]); % we choose a sligly higher number than the rank.
R = randn(l,m); % l x m
Y = R*A; % l x n
[q r Pr] = qr(Y); % q:l x l,    r: l x n,    Pr: n x n
r11 = r(1:k,1:k); %  k x k
r12 = r(1:k,k+1:end); % k x (n-k)
T = r11\r12; % k x (n-k)
P = [eye(k) T]*Pr.'; % P: k x n    
B = A*Pr(:,1:k);
s = Pr(:,1:k).'*(1:n).';