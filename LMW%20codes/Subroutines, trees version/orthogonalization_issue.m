% Script name:  orthogonalization_issue

% upload data
u_struct = load('beforeSVD');
u = u_struct.u;

% now the detailed space
u_detailed_struct = load('detailed');
u_detailed = u_detailed_struct.unit;

% u_detailed is orthogonal
norm(u_detailed'*u_detailed-eye(size(u_detailed,2)),'fro')

% and indeed, u \subset Sp(u_detailed)
norm(u-u_detailed*u_detailed'*u)

[q ,~ ,~] = svd(u,0);
% q is orthogonal
norm(q'*q-eye(size(q,2)),'fro')

% note that Sp(u) = Sp(q)
norm(u-q*q'*u)

% and, here is the punch line
norm(q-u_detailed*u_detailed'*q)

% Where is the numerical exactness disapear to???