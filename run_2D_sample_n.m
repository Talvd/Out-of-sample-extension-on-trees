%
% script to run all the samples of 2D extension in different algorithms and
% basis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

display_figure = 0;
% warning('off','all')
close all

new_angle = linspace(-pi+0.1,pi-0.1,150);
new_points =[cos(new_angle);sin(new_angle)];

%create the data for 2D samples
for N =[150 500 1000]

angle = linspace(-pi,pi,N);
points = [cos(angle);sin(angle)];

funcs = {};
%1
funcs{1}{1} = sin((points'*points(:,25) + 0.3).^(-1));
funcs{1}{2} = sin((new_points'*points(:,25) + 0.3).^(-1));
funcs{1}{3} = 'ocilla';

%2
f =  cos(2*points);
funcs{2}{1} = sqrt(sum(f.*f))';
new_f =  cos(2*new_points);
funcs{2}{2} = sqrt(sum(new_f.*new_f))';
funcs{2}{3} = 'smooth1';

%3
x0 = [1,0]';
shift = points-repmat(x0,1,size(points,2));
normshiftsqr = sum(shift.*shift);
funcs{3}{1} = exp(-normshiftsqr)';
new_shift = new_points-repmat(x0,1,size(new_points,2));
new_normshiftsqr = sum(new_shift.*new_shift);
funcs{3}{2} = exp(-new_normshiftsqr)';
funcs{3}{3} = 'smooth2';

%4
s = (points'*points(:,10))>0.1;
f = ones(size(points,2),1);
f(s) = .5;
funcs{4}{1} = f; 
new_s = (new_points'*points(:,10))>0.1;
new_f = ones(size(new_points,2),1);
new_f(new_s) = .5;
funcs{4}{2} = new_f; 
funcs{4}{3} = 'compactly_supported';

%5
r =[0.5    0.1;    0.8    0.3];
s = (points'*points(:,10))>0.1;
g = diag( points'*r*points);
f = zeros(size(points,2),1);
f(s) = g(s);
funcs{5}{1} = f; 
new_s = (new_points'*points(:,10))>0.1;
new_g = diag(new_points'*r*new_points);
new_f = zeros(size(new_points,2),1);
new_f(new_s) =new_g(new_s);
funcs{5}{2} = new_f;
funcs{5}{3}='compactly_supported_2';


funcs_num = size(funcs,2);

for i=1:funcs_num
    func = funcs{i};
    f = func{1};
    f_new_org = func{2};
    str = func{3};
    
    %LMW
    k_params = [1 10 15 N];
    
%     k_params =15;
    
    for k=k_params
        f_new =  run_LMW_extension(f, points, new_points, k, 2, [], 1);
        figure_title =['LMW k=' num2str(k) ' N=' num2str(N) ' function=' str];       
        figure_and_error_2D(points, new_points, f, f_new_org, f_new, figure_title, display_figure);
    end;
    
    %MSE
    
    data = [points new_points];
    d = squareform(pdist(data'));

    in_dists = d(1:N,1:N);
    out_dists = d(N+1:end,1:N);

    [AF S EF errors] = distancesMultiScale(f, in_dists, out_dists, 0.002 ,1);

    f_new_mse = EF(:,size(EF,2));
    figure_title =['MSE  N=' num2str(N) ' function=' str];   
    figure_and_error_2D(points, new_points, f, f_new_org, f_new_mse, figure_title, display_figure);
    
    %nystrom
    
    f_new =  run_nystrom_extension(f, points, new_points);
    figure_title =['nystrom  N=' num2str(N) ' function=' str];       
    figure_and_error_2D(points, new_points, f, f_new_org, f_new, figure_title, display_figure);
    
end
end

