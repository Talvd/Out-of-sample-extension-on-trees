%
% script to run all the samples of 2D extension in different algorithms and
% basis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

display_figure = 0;
warning('off','all')

%create the data for 2D samples
N = 150;

angle = linspace(-pi,pi,N);
points = [cos(angle);sin(angle)];

delta = abs((angle(1)-angle(2))*.5);
new_angle = linspace(-pi+delta,pi+delta,2*N-1);
new_points =[cos(new_angle);sin(new_angle)];

funcs = {@sin @step_sin @sin_multi};
funcs_num = size(funcs,2);

for i=1:funcs_num
    func = funcs{i};
    f = func(angle)';
    f_new_org = func(new_angle)';
    
    %LMW
    k_params = [1 10 15 N];
    
    for k=k_params
        f_new =  run_LMW_extension(f, points, new_points, k, 2);
        figure_title =['LMW k=' num2str(k) ' function=' char(func)];       
        figure_and_error_2D(points, new_points, f, f_new_org, f_new, figure_title, display_figure);
    end;
    
    %MSE
    
    data = [points new_points];
    d = squareform(pdist(data'));

    in_dists = d(1:N,1:N);
    out_dists = d(N+1:end,1:N);

    [AF S EF errors] = distancesMultiScale(f, in_dists, out_dists, 0.002 ,1);

    f_new_mse = EF(:,size(EF,2));
    figure_title =['MSE function=' char(func)];   
    figure_and_error_2D(points, new_points, f, f_new_org, f_new_mse, figure_title, display_figure);
    
    %nystrom
    
    f_new =  run_nystrom_extension(f, points, new_points);
    figure_title =['nystrom function=' char(func)];       
    figure_and_error_2D(points, new_points, f, f_new_org, f_new, figure_title, display_figure);
    
end

