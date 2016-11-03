function figure_and_error_2D(points,new_points,f, f_new_org, f_new, figure_title, display_figure)

% err_norm = norm(f_new_org - f_new);
rel_norm = norm(f_new_org - f_new)/norm(f_new_org);

err_str = num2str(rel_norm,'%.4f');
if (strcmp(err_str, '0.0000'))
    err_str = num2str(rel_norm,'%10.4e');
end

figure_title = [figure_title ' relative error ' err_str];

fprintf([figure_title '\n']);

if (display_figure)
    figure
    scatter3(points(1,:),points(2,:),f,'.b');
    hold
    scatter3(new_points(1,:),new_points(2,:),f_new,'.g');
    figure_title = strrep(figure_title, '_', '\_');
    title(figure_title)
end



end