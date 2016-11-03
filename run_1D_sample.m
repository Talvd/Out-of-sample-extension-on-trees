%
% script to run all the samples of 1D extension in different algorithms and
% basis
%
%---------------------------------------------------------------
% Nov 2016
%---------------------------------------------------------------

display_figure = 1;
more_figure=1;

run_ls = 0;
run_abcr = 0;
run_lmw = 1;
run_nystrom =0;
run_mse = 0;

warning('off','all')

%create the data for 1D samples
N = 64;
a = 0;
b = 2;
    
points = linspace(a,b,N);
delta = abs((points(1)-points(2))*.5);
new_points = linspace(a+delta,b,2*N-1);
    
% funcs = {@sine_one_over_x @step @sin @sin_15x @step_sin_1D};
funcs = {@sine_one_over_x @sin @sin_15x};

funcs_num = size(funcs,2);

for i=1:funcs_num
    func = funcs{i};
    f = func(points)';
    f_new_org = func(new_points)';
    
    % Least square
    if (run_ls)      
        f_new = [];
        for p=new_points
            v=quaddric_fit_extension(points,p,f,N,1);   
            f_new =[f_new; v(end)];               
        end 
        figure_title =['LS function=' char(func)];       
        figure_and_error_1D(points, new_points, f, f_new_org, f_new, figure_title, display_figure, more_figure);
    end
    
    %ABCR
    if (run_abcr)
        f_new =  run_ABCR_extension(f, points, new_points);
        figure_title =['ABCR function=' char(func)];       
        figure_and_error_1D(points, new_points, f, f_new_org, f_new, figure_title, display_figure, more_figure);
    end
    
    %LMW
    if (run_lmw)
%         k_params = [1 3 5 10 N];

        k_params = 10;

        for k=k_params
            f_new =  run_LMW_extension(f, points, new_points, k, 1);
            figure_title =['LMW k=' num2str(k) ' function=' char(func)];       
            figure_and_error_1D(points, new_points, f, f_new_org, f_new, figure_title, display_figure, more_figure);
        end;
    end
    
    %nystrom
    if (run_nystrom)
        f_new =  run_nystrom_extension(f, points, new_points);
        figure_title =['nystrom function=' char(func)];       
        figure_and_error_1D(points, new_points, f, f_new_org, f_new, figure_title, display_figure, more_figure);

        %nystrom fix epsilon  
        f_new =  run_nystrom_extension(f, points, new_points, (pi^2)/(2^7));
        figure_title =['nystrom fix eps function=' char(func)];       
        figure_and_error_1D(points, new_points, f, f_new_org, f_new, figure_title, display_figure, more_figure);
    end
    
    %MSE
    if (run_mse)
        data = [points new_points];
        d = squareform(pdist(data'));

        in_dists = d(1:N,1:N);
        out_dists = d(N+1:end,1:N);

        [AF S EF errors] = distancesMultiScale(f, in_dists, out_dists, 0.002 ,1);

        f_new_mse = EF(:,size(EF,2));
        figure_title =['MSE function=' char(func)];   
        figure_and_error_1D(points, new_points, f, f_new_org, f_new_mse, figure_title, display_figure, more_figure);
    end
    
   
end