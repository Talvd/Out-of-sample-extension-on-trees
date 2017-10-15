function [h]=figure_and_error_2D(points,new_points,f, f_new_org, f_new, figure_title, display_figure,h)

% err_norm = norm(f_new_org - f_new);

f_new_org(f_new_org==0) = 0.00000001;
rel_norm = norm(f_new_org - f_new)/norm(f_new_org);

err_str = num2str(rel_norm,'%.4f');
if (strcmp(err_str, '0.0000'))
    err_str = num2str(rel_norm,'%10.4e');
end

% figure_title = [figure_title ' relative error ' err_str];

fprintf([figure_title ' relative error ' err_str '\n']);

if (display_figure)
%     figure
%     scatter3(points(1,:),points(2,:),f,'ob','LineWidth',1.5);
%     hold
%     scatter3(new_points(1,:),new_points(2,:),f_new,150,(f_new_org - f_new)./f_new_org,'Marker','.');
%     figure_title = strrep(figure_title, '_', '\_');
%     title(figure_title)
 
    rel_err = abs((f_new_org - f_new)./f_new_org);
    if nargin<8
        figure
        plot(rel_err);
    else
        figure(h) 
        [~,~,~,text]=legend;
        if (isempty(text))
            p = plot(rel_err);
            hold all
        elseif (size(text,2) == 1)
                p = plot(rel_err, ':');
%             set(p,'MarkerSize', 10)
        else
            p = plot(rel_err, '-.');
        end         
        
        set(p, 'LineWidth',1.5);
        ind =  strfind(figure_title, 'function');
        fig_legend = figure_title(1:ind-1);        
        legend(text,fig_legend)            
    end         
end