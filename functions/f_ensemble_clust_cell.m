function clust_out = f_ensemble_clust_cell(coeffs, scores, raster_norm, params)

ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'
shuffle_method = f_get_param(params, 'shuffle_method', 'circ_shift');
corr_cell_thresh_percent = f_get_param(params, 'corr_cell_thresh_percent', 95);

hcluster_method = f_get_param(params, 'hcluster_method', 'ward');
hcluster_distance_metric = f_get_param(params, 'hcluster_distance_metric', 'cosine');

plot_stuff = f_get_param(params, 'plot_stuff', 0);

[num_cells, num_clust] = size(coeffs);

% if ~strcmpi(params.hcluster_distance_metric, 'cosine')
%     % cosine will not pull out the no activity cluster
%     num_clust = num_clust + 1;
% end 

%% get thresh for correlated cells
dist_d = f_pdist_YS(raster_norm, hcluster_distance_metric);
dist_d2 = dist_d + diag(ones(num_cells,1))*mean(squareform(dist_d));
pwcorr_d = mean(1-dist_d2,2);

num_reps = 100;
pwcorr_s_means = zeros(num_cells,num_reps);
for n_rep = 1:num_reps
    raster_norm_shuff = f_shuffle_data(raster_norm, shuffle_method);
    dist_s = f_pdist_YS(raster_norm_shuff, hcluster_distance_metric);
    dist_s2 = dist_s + diag(ones(num_cells,1))*mean(squareform(dist_s));
    pwcorr_s_means(:,n_rep) = mean(1-dist_s2);
end

dist_s_thresh_up = prctile(pwcorr_s_means, corr_cell_thresh_percent,2);
dist_s_thresh_down = prctile(pwcorr_s_means, 100-corr_cell_thresh_percent,2);

if plot_stuff
    figure; hold on; axis tight;
    plot(pwcorr_d);
    plot(dist_s_thresh_up, '--r');
    plot(dist_s_thresh_down, '--g');
    legend('mean pw corr', [num2str(corr_cell_thresh_percent) '% shuff thresh'], [num2str(100-corr_cell_thresh_percent) '% shuff thresh'])
    xlabel('cells'); ylabel('mean pw corr');
    title('correlated cells selection');
end

corr_cells = pwcorr_d>dist_s_thresh_up;

%%

X = coeffs(corr_cells,:);

%% cluster cells
if strcmpi(cluster_method, 'hclust')
    %% cluster with hclust
    params2.method = hcluster_method; % ward(inner square), average, single(shortest)
    params2.distance_metric = 'euclidean'; %hcluster_distance_metric; % none, euclidean, squaredeuclidean, cosine, hammilarity, rbf
    params2.plot_dist_mat = plot_stuff;
    params2.plot_clusters = plot_stuff;
    params2.num_clust = num_clust;
    params2.XY_label = 'Cells';
    params2.clim = [0 1];
    clust_out_cell = f_hcluster_wrap(X, params2);
    %gscatter(X(:,1),X(:,2),hclust_out.clust_ident);
    
elseif strcmpi(cluster_method, 'gmm')
    %% clust with gmm, can find best regularizer
    optimize_reg_val = 0;
    if optimize_reg_val
        num_vals = 50;
        num_reps = 10;
        mean_acc = zeros(num_vals,num_reps);
        %rg_list = linspace(0.0001, 1, num_vals);
        rg_list = logspace(-4, -2, num_vals);
        for n_rg = 1:numel(rg_list)
            for n_rep = 1:num_reps
                params3.metric = 'cosine'; % cosine squaredeuclidean
                params3.RegularizationValue = rg_list(n_rg);
                params3.num_clust = num_ens+1;
                gmmclust_out = f_gmmcluster_trial(X, params3);
                %gscatter(X(:,1),X(:,2),gmmclust_out.clust_ident);
                eval_gmm = f_evaluate_ens_result(gmmclust_out.clust_ident, ens_list_gt, 0);
                mean_acc(n_rg, n_rep) = mean(eval_gmm.accuracy);
            end
        end
        f_rg = figure;
        shadedErrorBar(log10(rg_list), mean(mean_acc,2), std(mean_acc, [], 2)/sqrt(num_reps-1));
        xlabel('log10(rg val)')
        title('Click on optimal regression value');
        [x,~] = ginput(1);
        close(f_rg);
        params3.RegularizationValue = 10^x;
    else
        params3.RegularizationValue = 0.006;
    end
    params3.metric = 'sqEeuclidean'; % cosine squaredeuclidean
    params3.num_clust = num_clust;
    clust_out_cell = f_gmmcluster_trial(X, params3);
    [~, clust_out_cell.dend_order] = sort(clust_out_cell.clust_ident);
end

%% put back all uncorrelated cells

clust_out_cell_full = clust_out_cell;
clust_out_cell_full.corr_cells = corr_cells;
clust_out_cell_full.clust_out_cell_corr = clust_out_cell;

clust_ident_full = zeros(num_cells,1);
clust_ident_full(corr_cells) = clust_out_cell.clust_ident;
clust_out_cell_full.clust_ident = clust_ident_full;

Xn = coeffs(~corr_cells,:);
clust_out_cell_uncorr = f_hcluster_wrap(Xn, params2);

cell_list = (1:num_cells)';
cell_list_X = cell_list(corr_cells);
cell_list_Xn = cell_list(~corr_cells);
den_order_full = [cell_list_X(clust_out_cell.dend_order); cell_list_Xn(clust_out_cell_uncorr.dend_order)];
clust_out_cell_full.dend_order = den_order_full;

clust_out_cell_full2 = f_get_clust_params(coeffs, clust_out_cell_full);

% dist1 = f_pdist_YS(X(clust_out_cell.dend_order,:), 'cosine');
% figure; imagesc(1-dist1)
% 
% dist1 = f_pdist_YS(Xn(clust_out_cell_uncorr.dend_order,:), 'cosine');
% figure; imagesc(1-dist1)
% 
% dist1 = f_pdist_YS(coeffs(den_order_full,:));
% figure; imagesc(1-dist1)

num_frames = size(scores,2);
% get score projections
ens_scores = zeros(num_clust,num_frames);
ens_fames = cell(num_clust,1);
for n_cl = 1:num_clust
    score1 = coeffs(clust_out_cell_full2.clust_ident == n_cl,:)*scores;
    if numel(clust_out_cell_full2.ens_list{n_cl})>1
        score1 = mean(score1);
    end
    ens_scores(n_cl,:) = score1;
    ens_fames{n_cl} = find(score1>(3*std(score1)));
end
clust_out_cell_full2.ens_scores = ens_scores;

clust_out.cells = clust_out_cell_full2;
clust_out.trials.ens_list = ens_fames;

% for n_cl = 1:num_clust
%     ens_raster = raster_norm(clust_out_cell_full2.clust_ident == n_cl,:);
%     figure;
%     subplot(2,1,1);
%     imagesc(ens_raster);
%     title(['clust ' num2str(n_cl)])
%     subplot(2,1,2); hold on;
%     plot(ens_scores(n_cl,:)); axis tight;
%     plot(ones(num_frames,1)*3*std(ens_scores(n_cl,:)), '--r')
% end


end