function s = ColByColRandomInterpolativeDecomposition(rows_number,cols_number, k, generate_mat_column,data)

% Columnwise random  interpolative decomposition algorithm%

% INPUT:
% rows_number - number of rows of the decomposed matrix A
% cols_number - number of columns of the decomposed matrix A
% k - the number of columns to return. The numerical rank
% generate_mat_columns(data,j) - function handle. that generates the j-th column of A, using the parameter data

% OUTOUT:
% s - the chosen columns of A.

% NOTE: 1. If we multiply B and P we get an approximation to A.  A ~ B*P 
%       2. If only the selected indexes are needed (s), comment out lines
%       below and use [s, ~, ~] =
%       ColByColRandomInterpolativeDecomposition(...)



m = rows_number;
n = cols_number;

l = min([k+3 m n]); % we choose a slightly higher number than the rank.
R = randn(l,m); % l x m

Y = zeros(l, n);

for ii=1:n
    col = generate_mat_column(data,ii);
    Y(:,ii) = R*col;
end

[~, ~, Pr] = qr(Y); % q:l x l,    r: l x n,    Pr: n x n
s = Pr(:,1:k).'*(1:n).';

