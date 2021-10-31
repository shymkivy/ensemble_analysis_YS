function ens_out = f_ensemble_analysis_YS_raster(firing_rate, params)
% input either 3D tials data (Cell x Time X Trial)
%           or 2D trial data (Cell x Trial)
%% parameters
if ~exist('params', 'var') || ~isstruct(params)
    params = struct;
end

normalize1 = f_get_param(params, 'normalize', 'norm_mean_std'); % 'norm_mean_std', 'norm_mean' 'none'
shuffle_method = f_get_param(params, 'shuffle_method', 'circ_shift');     % 'circ_shift' or 'scramble'
total_dim_thresh = f_get_param(params, 'total_dim_thresh', .7);
dim_est_num_reps = f_get_param(params, 'dim_est_num_reps', 50);
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf'); % 'PCA', 'AV', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
ensamble_extraction = f_get_param(params, 'ensamble_extraction', 'thresh'); % clust 'thresh'
plot_stuff = f_get_param(params, 'plot_stuff', 0);

num_comps = f_get_param(params, 'num_comp');

if ~strcmpi(ensamble_method, 'nmf') && strcmpi(ensamble_extraction, 'thresh')
    fprintf('Thresh detection is only for nmf dummy...\n');
    ensamble_extraction = 'clust';
end

fprintf('Detecting ensembles with %s and %s...',ensamble_method, ensamble_extraction);

%%
ndims1 = ndims(firing_rate);
if ndims1 == 3
    num_cells = size(firing_rate,1);
    firing_rate = reshape(firing_rate, num_cells,[]);
end

% active_cells = sum(firing_rate,2) ~= 0;
% firing_rate(~active_cells,:) = [];
% ens_out.active_cells = active_cells;

firing_rate_norm = f_normalize(firing_rate, normalize1);

[num_cells, num_T] = size(firing_rate_norm);

%% dim reduction with SVD to calulate components number

% [U,S,V] = svd(firing_rate_norm);
% sing_val_sq = diag(S'*S);
% d_explained = sing_val_sq/sum(sing_val_sq)*100;

[d_coeff,d_score,~,~,d_explained,d_mu] = pca(firing_rate_norm');

dimensionality_total_norm = sum(cumsum(d_explained)<(total_dim_thresh*100));
%figure; plot(d_explained)

%% repeat with not norm
[~,S2,~] = svd(firing_rate);
sing_val_sq2 = diag(S2'*S2);
d_explained2 = sing_val_sq2/sum(sing_val_sq2)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained2)<(total_dim_thresh*100));


ens_out.dimensionality_total = dimensionality_total;
ens_out.dimensionality_total_norm = dimensionality_total_norm;
ens_out.dimensionality_first_comp_size = d_explained2(1);
ens_out.d_explained = d_explained(1:num_comps);

%% shuff and PCA
if isempty(num_comps)
    fprintf('Estimating dimensionality...\n');
    dim_est_num_reps = 50;
    max_lamb_shuff = zeros(dim_est_num_reps,1);
    dim_total_shuff = zeros(dim_est_num_reps,1);
    for n_rep = 1:dim_est_num_reps
        firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
    %     [~,s_S,~] = svd(firing_rate_shuff);
    %     s_sing_val_sq = diag(s_S'*s_S);
    %     s_explained = s_sing_val_sq/sum(s_sing_val_sq)*100;
        [~,~,~,~,s_explained,~] = pca(firing_rate_shuff');

        dim_total_shuff(n_rep) = sum(cumsum(s_explained)<(total_dim_thresh*100));
        max_lamb_shuff(n_rep) = max(s_explained);
    end
    dimensionality_total_norm_shuff = mean(dim_total_shuff);
    % eigenvalues below lower bound plus above upper should
    % theoretically equal total number of neurons in all ensembles
    %dimensionality_corr = sum(d_explained>prctile(max_lamb_shuff, corr_comp_thresh*100));

    dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'));
    num_comps = ceil(dimensionality_corr);
    
    ens_out.dimensionality_total_norm_shuff = dimensionality_total_norm_shuff;
    ens_out.dimensionality_corr = dimensionality_corr;
    fprintf('; estimated %d dimensions...\n',num_comps);
else
    fprintf('; using %d components\n', num_comps);
end
ens_out.num_comps = num_comps;

%%
%SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
%SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);

%firing_rate_LR = U(:,1:num_comps)*S(1:num_comps,1:num_comps)*V(:,1:num_comps)';
%SI_firing_rate_LR = similarity_index(firing_rate_LR, firing_rate_LR);

n_comp = 1:num_comps;
firing_rate_LR = (d_coeff(:,n_comp)*d_score(:,n_comp)'+d_mu')';



%% sort cells and trials
hc_params.method = params.hcluster_method;
hc_params.distance_metric = params.hcluster_distance_metric;
hc_params.plot_dist_mat = plot_stuff;
hc_params.plot_clusters = plot_stuff;
hc_params.num_clust = num_comps+1;
hc_params.title_tag = 'Coeffs (cells)';
hclust_out_cell = f_hcluster_wrap(d_coeff(:,n_comp), hc_params);
ord_cell = hclust_out_cell.dend_order;
ens_out.ord_cell = ord_cell;

sort_tr = 0;
if sort_tr
    d_score_norm = d_score(:,n_comp)./vecnorm(d_score(:,n_comp));
    hc_params.title_tag = 'Scores (trials)';
    hclust_out_tr = f_hcluster_wrap(d_score_norm, hc_params);
    ord_tr = hclust_out_tr.dend_order;
    ens_out.ord_tr = ord_tr;
end
%% real data 
if num_comps > 0
    fprintf('Dim reduction with %s, %d comps...\n', ensamble_method, num_comps);
    num_ens_comps = num_comps;
    firing_rate_ensemb = firing_rate_norm;
    
    [dred_factors1, ~] = f_dred_train2(firing_rate_ensemb, num_ens_comps, ensamble_method, 0);
    [coeffs, scores] = f_dred_get_coeffs(dred_factors1);

    ens_out.ensamble_method = ensamble_method;
    ens_out.num_ens_comps = num_ens_comps;
    ens_out.dred_factors = dred_factors1;
    ens_out.coeffs = coeffs;
    ens_out.scores = scores;
    %%
    %disp('Extraction...')
    if strcmpi(ensamble_extraction, 'clust')
        ens_out1 = f_ensemble_clust_cell(coeffs, scores, firing_rate_norm, params);
        %ens_out1 = f_ensemble_extract_clust(coeffs, scores, num_ens_comps, firing_rate_norm, params);
    elseif strcmpi(ensamble_extraction, 'thresh')
        [thresh_coeffs, thresh_scores] = f_ens_get_thresh(firing_rate_ensemb, coeffs, scores, num_ens_comps, params);
        ens_out1 = f_ensemble_apply_thresh(coeffs, scores, thresh_coeffs, thresh_scores, num_ens_comps);
        ens_out1.cells.ens_scores = ens_out.scores;
    end
    ens_out.cells = ens_out1.cells;
    ens_out.trials = ens_out1.trials;
else
    ens_out.cells.clust_label = 0;
    ens_out.cells.residual_list = {(1:num_cells)'};
    ens_out.cells.clust_ident = zeros(num_cells,1);
    ens_out.trials.clust_label = 0;
    ens_out.trials.residual_list = {(1:num_T)'};
    ens_out.trials.clust_ident = zeros(num_trials,1);

end
ens_out.cells.dend_order = ord_cell;
if sort_tr
    ens_out.trials.dend_order = ord_tr;
end

%%
if plot_stuff
    f_plot_raster_mean(firing_rate, 1)
    title('raster pre')

    f_plot_raster_mean(firing_rate(ord_cell,:),1)
    title('raster cell sort')
    
    n_comp = 1:num_ens_comps;
    firing_rate_LR2 = coeffs(:,n_comp)*scores(n_comp,:);
    f_plot_raster_mean(firing_rate_LR2(ord_cell,:), 1)
    title(['raster ' ensamble_method ' LR2 cell sort'])
        
    if sort_tr
        f_plot_raster_mean(firing_rate_norm(ord_cell,ord_tr), 1)
        title('raster cell trial sort')

        f_plot_raster_mean(firing_rate_LR(ord_cell,ord_tr), 1)
        title('raster SVD LR cell trial sort')

        f_plot_raster_mean(firing_rate_LR2(ord_cell,ord_tr), 1)
        title(['raster ' ensamble_method ' LR2 cell trial sort'])
    end
end


end
