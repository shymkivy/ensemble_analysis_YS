function num_clusters = f_estimate_clust_num(data, params)

method = f_get_param(params, 'method');
metric = f_get_param(params, 'metric');
num_shuff = f_get_param(params, 'num_shuff', 10);
clust_range = f_get_param(params, 'clust_range', 1:8);
plot_stuff = f_get_param(params, 'plot_stuff', 0);
plot_stuff_full = f_get_param(params, 'plot_stuff_full', 0);
manual_cluster_input = f_get_param(params, 'manual_cluster_input', 0);
title_tag = f_get_param(params, 'title_tag');

[num_groups, num_features] = size(data);

negative_data = sum(data(:)<0);

% first estimate within clust variance for shuff data
mean_clust_var_shuff = zeros(numel(clust_range),num_shuff);
for n_shuff = 1:num_shuff
    data_shuff = data;
    % shiffle within each axis
    vec1 = zeros(num_features,1);
    for n_vec = 1:num_groups
        % generate random angle in positive n dim data space and rotate
        % data vec to that angle
        angs = rand(num_features-1,1)*pi/(2-negative_data);
        r = 1;
        curr_r = r;
        for n_dim = 1:(num_features-1)
            vec1(n_dim) = curr_r * sin(angs(n_dim));
            curr_r = curr_r * cos(angs(n_dim));
        end
        vec1(end) = curr_r;        
        %[vec1(3),vec1(2),vec1(1)] = sph2cart(angs(2),angs(1), 1);
        data_shuff(n_vec,:) = norm(data_shuff(n_vec,:))*vec1;
    end
%     for n_col = 1:num_features
%         data_shuff(:,n_col) = data_shuff(randperm(num_groups),n_col);
%     end
    % get var
    for num_clust1 = clust_range
        [~, clust_ident_shuff, ~] = f_hcluster(data_shuff, method, num_clust1);
        mean_clust_var_shuff(num_clust1,n_shuff) = if_get_cluster_dist_var(data_shuff, clust_ident_shuff, metric, num_clust1);
        if n_shuff == 1
            if plot_stuff_full
                f_plot_comp_scatter(data_shuff, clust_ident_shuff);
                title([title_tag '; shuff']);
            end
        end
    end
    
end
% and real data
mean_clust_var = zeros(numel(clust_range),1);
for num_clust1 = clust_range
    [~, clust_ident, ~] = f_hcluster(data, method, num_clust1);
    mean_clust_var(num_clust1) = if_get_cluster_dist_var(data, clust_ident, metric, num_clust1);
    if plot_stuff_full
        f_plot_comp_scatter(data, clust_ident);
        title([title_tag '; data']);
    end
end

if plot_stuff
    shuff1 = mean(mean_clust_var_shuff,2);
    shuff1_sem = std(mean_clust_var_shuff,[],2)/sqrt(max(num_shuff-1,1));
    
    shuff2 = shuff1/shuff1(1);
    shuff2_sem = std([zeros(1,numel(clust_range)); diff(mean_clust_var_shuff/shuff1(1))],[],2)/sqrt(max(num_shuff-1,1));
    data1 = mean_clust_var/mean_clust_var(1);
    
    figure; 
    subplot(2,1,1); hold on;
    plot(clust_range, mean_clust_var, 'LineWidth', 2)
    errorbar(clust_range, shuff1, shuff1_sem, 'LineWidth', 2, 'Color', [.4 .4 .4])
    title([title_tag '; Total error']);
    legend('Data', 'Shuff');
    
    subplot(2,1,2); hold on;
    plot(clust_range, [0; diff(data1)], 'LineWidth', 2)
    errorbar(clust_range, [0; diff(shuff2)], shuff2_sem, 'LineWidth', 2, 'Color', [.4 .4 .4])
    title('Fraction error change with new components');
    legend('Data', 'Shuff')
end

if manual_cluster_input
    num_clusters = input('Rafa: How many clusters are there?');
else
    num_clusters = NaN;
end
end


function cluster_var = if_get_cluster_dist_var(data, clust_ident, metric, num_clust)

dist_cell = cell(num_clust,1);
for n_clust = 1:num_clust
    temp_data = data(clust_ident == n_clust,:);
    mean_vec = mean(temp_data);
    dist_cell{n_clust} = f_pdist2_YS(temp_data, mean_vec, metric);
end
dist_all = cat(1,dist_cell{:});
cluster_var = mean(dist_all);

end