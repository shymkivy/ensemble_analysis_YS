function hclust_out = f_hcluster_wrap(X, params)
% input X (samples x features)
% clusters into types of samples

num_clust = f_get_param(params, 'num_clust');
estimate_clust_num = f_get_param(params, 'estimate_clust_num', 0);
method = f_get_param(params, 'method', 'cosine');
metric = f_get_param(params, 'metric', 'cosine');
sp = f_get_param(params, 'subplot_ptr');
plot_dist_mat = f_get_param(params, 'plot_dist_mat', 1);
plot_clusters = f_get_param(params, 'plot_clusters', 1);
XY_label = f_get_param(params, 'XY_label');
title_tag = f_get_param(params, 'title_tag');
sample_types = f_get_param(params, 'sample_types');
sample_types_colors = f_get_param(params, 'sample_types_colors');

if isempty(num_clust)
    num_clust = 1;
end

% warning is because some trial bins are nearly zero
[dend_order, clust_ident, Z] = f_hcluster(X, method, num_clust);

%f_plot_comp_scatter(trial_peaks, clust_ident)

%figure; imagesc(color_seq_temporal)

dist1 = f_pdist_YS(X(dend_order,:), metric);

if estimate_clust_num
    est_params.num_shuff = 5;
    est_params.clust_range = 1:5;
    est_params.method = method;
    est_params.metric = metric;
    est_params.manual_cluster_input = 0;
    est_params.plot_stuff = 1;
    est_params.title_tag = title_tag;
    f_estimate_clust_num(X, est_params);
end

hclust_out.dist = dist1;
hclust_out.dend_order = dend_order;
hclust_out.clust_ident = clust_ident;
hclust_out.Z = Z;

if plot_dist_mat
    image_Z = 1-dist1;
    if isempty(sp)
        figure;
        sp = gca;
    else
        subplot(sp);
    end
    imagesc(image_Z);
    %axis image;
    title(sprintf('%s %s hclust, %s dist',  title_tag, method, metric));
    caxis([0 1]);
    axis tight;
    axis equal;
    if ~isempty(XY_label)
        xlabel(XY_label);
        ylabel(XY_label);
    end
    %colorbar;
    clim1 = caxis;
    sp.YDir = 'reverse';

    %% add variable type indicator

    if ~isempty(sample_types)
        f_plot_trial_indicator(sample_types, dend_order, 1, numel(sample_types), sample_types_colors);
    end
    %imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
    %imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

    %% plot clusters
    if plot_clusters
        ord1 = clust_ident(dend_order);
        for num_clust1 = 1:num_clust
            temp_list = find(ord1 == num_clust1);
            subplot(sp); hold on;
            rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
        end
    end
    hclust_out.clim = clim1;
end

end
