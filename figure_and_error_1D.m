function [h] = figure_and_error_1D(points,new_points,f, f_new_org, f_new, in_figure_title, display_figure, more_figures, h)

% err_norm = norm(f_new_org - f_new);
rel_norm = norm(f_new_org - f_new)/norm(f_new_org);

err_str = num2str(rel_norm,'%.4f');
if (strcmp(err_str, '0.0000'))
    err_str = num2str(rel_norm,'%10.4e');
end

figure_title = [in_figure_title ' relative error ' err_str];

fprintf([figure_title '\n']);

if (display_figure)
    figure
    scatter(new_points,f_new,'.g');
    hold
    scatter(points,f,'.b');
    figure_title = strrep(figure_title, '_', '\_');
    title(figure_title)
    if(more_figures)
%         err = abs(f_new_org - f_new);
%         figure
%         plot(new_points,err);
%         legend('ERROR');
        rel_err = abs((f_new_org - f_new)./f_new_org);
        if nargin<9
            figure
            plot(new_points,rel_err);
        else
            figure(h)
            hold all
            p = plot(new_points,rel_err);
            set(p, 'LineWidth',2);
            ind =  strfind(in_figure_title, 'function');
            fig_legend = in_figure_title(1:ind-1);
            [~,~,~,text]=legend;
            legend(text,fig_legend)            
        end               
        title('Relative ERROR');
        figure
        d = diff(f')./diff(points);
        plot(points(1:(end-1)),d);
        legend('Derivative');
    end
end

end
