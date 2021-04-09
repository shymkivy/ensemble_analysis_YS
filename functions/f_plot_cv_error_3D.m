function f_plot_cv_error_3D(data, data_s, x_var, y_var, z_var)

data_x = [data.(x_var)];
data_y = [data.(y_var)];
data_z = [data.(z_var)];
data_reps = [data.n_rep] ;

data_all = [data_x(:), data_y(:), data_reps(:), data_z(:)];

[~, sort_ind1] = sort(data_all(:,1), 'ascend');
[~, sort_ind2] = sort(data_all(sort_ind1,2), 'ascend');
[~, sort_ind3] = sort(data_all(sort_ind1(sort_ind2),3), 'ascend');

data_all2 = data_all(sort_ind1(sort_ind2(sort_ind3)),:);

numel_x = numel(unique(data_x));
numel_y = numel(unique(data_y));
numel_rep = numel(unique(data_reps));

grid_x = reshape(data_all2(:,1), [numel_x, numel_y, numel_rep]);
grid_y = reshape(data_all2(:,2), [numel_x, numel_y, numel_rep]);
grid_z = reshape(data_all2(:,4), [numel_x, numel_y, numel_rep]);

figure; hold on;
if numel_x == 1
    p1 = plot(squeeze(grid_y(:,:,1)), mean(grid_z,3), 'b', 'LineWidth', 2);
    if numel_rep > 1
        errorbar(squeeze(grid_y(:,:,1)), mean(grid_z,3), std(grid_z,[],3), 'b', 'LineWidth', 2)
    end
    xlabel(y_var, 'interpreter', 'none');
    ylabel(z_var, 'interpreter', 'none');
elseif numel_y == 1
    p1 = plot(squeeze(grid_x(:,:,1)), mean(grid_z,3), 'b', 'LineWidth', 2);
    if numel_rep > 1
        errorbar(squeeze(grid_x(:,:,1)), mean(grid_z,3), std(grid_z,[],3), 'b', 'LineWidth', 2)
    end
    xlabel(x_var, 'interpreter', 'none');
    ylabel(z_var, 'interpreter', 'none');
else
    p1 = surf(squeeze(grid_x(:,:,1)), squeeze(grid_y(:,:,1)), mean(grid_z,3));
    xlabel(x_var, 'interpreter', 'none');
    ylabel(y_var, 'interpreter', 'none');
    zlabel(z_var, 'interpreter', 'none');
end
title(z_var, 'interpreter', 'none')

if ~isempty(data_s)
    data_x = [data_s.(x_var)];
    data_y = [data_s.(y_var)];
    data_z = [data_s.(z_var)];
    data_reps = [data_s.n_rep] ;

    data_all = [data_x(:), data_y(:), data_reps(:), data_z(:)];

    [~, sort_ind1] = sort(data_all(:,1), 'ascend');
    [~, sort_ind2] = sort(data_all(sort_ind1,2), 'ascend');
    [~, sort_ind3] = sort(data_all(sort_ind1(sort_ind2),3), 'ascend');

    data_all2 = data_all(sort_ind1(sort_ind2(sort_ind3)),:);

    numel_x = numel(unique(data_x));
    numel_y = numel(unique(data_y));
    numel_rep = numel(unique(data_reps));

    grid_x = reshape(data_all2(:,1), [numel_x, numel_y, numel_rep]);
    grid_y = reshape(data_all2(:,2), [numel_x, numel_y, numel_rep]);
    grid_z = reshape(data_all2(:,4), [numel_x, numel_y, numel_rep]);

    if numel_x == 1
        p2 = plot(squeeze(grid_y(:,:,1)), mean(grid_z,3), 'k', 'LineWidth', 2);
        if numel_rep > 1
            errorbar(squeeze(grid_y(:,:,1)), mean(grid_z,3), std(grid_z,[],3), 'k', 'LineWidth', 2)
        end
        xlabel(y_var, 'interpreter', 'none');
        ylabel(z_var, 'interpreter', 'none');
    elseif numel_y == 1
        p2 = plot(squeeze(grid_x(:,:,1)), mean(grid_z,3), 'k', 'LineWidth', 2);
        if numel_rep > 1
            errorbar(squeeze(grid_x(:,:,1)), mean(grid_z,3), std(grid_z,[],3), 'k', 'LineWidth', 2)
        end
        xlabel(x_var, 'interpreter', 'none');
        ylabel(z_var, 'interpreter', 'none');
    else
        p2 = surf(squeeze(grid_x(:,:,1)), squeeze(grid_y(:,:,1)), mean(grid_z,3));
        xlabel(x_var, 'interpreter', 'none');
        ylabel(y_var, 'interpreter', 'none');
        zlabel(z_var, 'interpreter', 'none');
    end
    title(z_var, 'interpreter', 'none') 
end

legend([p1 p2], {'data', 'shuff'})

end