% restriction_script

f0 = ones(20,2);
f0(:,2) = 1:20;

f0 = f0*diag(sum(f0.*f0).^(-.5));

f1 = zeros(20,2);
f2 = zeros(20,2);
f1(1:10,:) = f0(1:10,:);
f2(11:20,:) = f0(11:20,:);

f1 = f1*diag(sum(f1.*f1).^(-.5));
f2 = f2*diag(sum(f2.*f2).^(-.5));

[q ,~ ,~] = svd([f1,f2],0);

norm(f2-q*q'*f2)
norm(f1-q*q'*f1)
norm(f0-q*q'*f0)
