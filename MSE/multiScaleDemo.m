% This is a short script to explain how to use distancesMultiScale()
%% generate input
% n = 100;
% x = linspace(0,2*pi,n)';
% f = sin(x.*x);

n = 64;
% a = 0;
% b = 2;
% x = linspace(a,b,n)';
% delta = abs((x(1)-x(2))*.5);
% y = linspace(a+delta,b,2*n-1)';

func = @sin_15x; %@step;sin_15x;sine_one_over_x;sin

f = func(x);

% y = linspace(0,2*pi,2*n)';
data = [x ; y];
d = squareform(pdist(data));

in_dists = d(1:n,1:n);
out_dists = d(n+1:end,1:n);

%% calculate the approximation
[AF S EF errors] = distancesMultiScale(f, in_dists, out_dists, 0.002 ,1);

%% draw results
figure
hold on;
scatter(y,EF(:,1),'g')
scatter(y,EF(:,size(EF,2)),'b')
plot(y,EF(:,size(EF,2)))
scatter(x,f,'r')
legend('First iteration Approximation','Last Iteration Approximation',...
'Last Iteration Approximation (Line)','Original Samples','Location','SouthWest');

f_new_org = func(y);
f_new = EF(:,size(EF,2));
err_norm = norm(f_new_org - f_new)
rel_norm = norm(f_new_org - f_new)/norm(f_new_org)

figure
scatter(y,f_new,'.g');
hold
h = scatter(x,f,'.b');
