function ens_out = f_ensemble_extract_clust(coeffs, scores, num_ens, params)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'
cluster_method_cell = f_get_param(params, 'cluster_method_cell', 'hclust');
plot_stuff = f_get_param(params, 'plot_stuff', 0);

num_comps = size(scores,1);
num_clust = num_ens + 1;

%% cluster trials from scores

X = scores';

if strcmpi(cluster_method, 'hclust')
    %% cluster with hclust
    params2.method = 'ward'; % cosine, ward
    params2.metric = 'squaredeuclidean'; % cosine squaredeuclidean
    params2.plot_dist_mat = plot_stuff;
    params2.plot_clusters = plot_stuff;
    params2.num_clust = num_clust;
    params2.XY_label = 'Trials';
    clust_out_tr = f_hcluster_wrap(X, params2);
    %gscatter(X(:,1),X(:,2),hclust_out.clust_ident);
    
elseif strcmpi(cluster_method, 'gmm')
    %% clust with gmm, can find best regularizer
    optimize_reg_val = 0;
    if optimize_reg_val
        num_vals = 100;
        num_reps = 20;
        mean_acc = zeros(num_vals,num_reps);
        %rg_list = linspace(0.0001, 1, num_vals);
        rg_list = logspace(-5, 0, num_vals);
        for n_rg = 1:numel(rg_list)
            for n_rep = 1:num_reps
                params3.metric = 'cosine'; % cosine squaredeuclidean
                params3.RegularizationValue = rg_list(n_rg);
                params3.num_clust = num_comps+1;
                gmmclust_out = f_gmmcluster_trial(X, params3);
                %gscatter(X(:,1),X(:,2),gmmclust_out.clust_ident);
                eval_gmm = f_evaluate_ens_result(gmmclust_out.clust_ident, trial_list_gt, 0);
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
        params3.RegularizationValue = 0.017;
    end
    params3.metric = 'cosine'; % cosine squaredeuclidean
    params3.num_clust = num_clust;
    clust_out_tr = f_gmmcluster_trial(X, params3);
    [~, clust_out_tr.dend_order] = sort(clust_out_tr.clust_ident);
end

clust_params_tr = if_get_clust_params(X, clust_out_tr);


%% get the ensemble cells from trials
  
% proj_score = scores'/(scores*scores');
% %proj_clust = clust_cent'/(clust_cent*clust_cent');
% 
% %cell_proj = coeffs*proj_clust;
% 
% %coeffs2 = firing_rate_norm*proj
% clust_cent = ens_out.clust_centers_tr;
% 
% num_reps = 20;
% proj_mags = zeros(num_cells, num_reps, num_clust);
% for n_rep = 1:num_reps
%     firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
%     coeffs_s = firing_rate_shuff*proj_score;
%     for n_cl = 1:num_clust
%         proj_mags(:,n_rep,n_cl) = coeffs_s*clust_cent(n_cl,:)'/norm(clust_cent(n_cl,:));
%     end
% end
% proj_mags = reshape(proj_mags,[],num_clust);
% 
% ens_out.ens_cells = cell(num_clust,1);
% ens_out.ens_cell_z_prob = cell(num_clust,1);
% z_thresh = 2;
% for n_cl = 2:num_clust
%     clust_cell_proj = coeffs*clust_cent(n_cl,:)'/norm(clust_cent(n_cl,:));
%     
%     mean1 = mean(proj_mags(:,n_cl));  
%     z_fac = prctile(proj_mags(:,n_cl)- mean1,95)/2;
%     
%     clust_cell_proj_z = (clust_cell_proj-mean1)/z_fac;
%     [cell_z, cell_seq] = sort(clust_cell_proj_z, 'descend');
%     
%     ens_out.ens_cells{n_cl} = cell_seq(cell_z>z_thresh);
%     ens_out.ens_cell_z_prob{n_cl} = cell_z(cell_z>z_thresh);
% end

%% try cluster coeffs for cells

X = coeffs;
%X = firing_rate_LR;

%% cluster cells
if strcmpi(cluster_method_cell, 'hclust')
    %% cluster with hclust
    params2.method = 'ward'; % cosine, ward
    params2.metric = 'squaredeuclidean'; % cosine squaredeuclidean
    params2.plot_dist_mat = plot_stuff;
    params2.plot_clusters = plot_stuff;
    params2.num_clust = num_clust;
    params2.XY_label = 'Cells';
    clust_out_cell = f_hcluster_wrap(X, params2);
    %gscatter(X(:,1),X(:,2),hclust_out.clust_ident);
    
elseif strcmpi(cluster_method_cell, 'gmm')
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
                params3.num_clust = num_comps+1;
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

clust_params_cell = if_get_clust_params(X, clust_out_cell);

%% align to trial clusts


[cells_clust_al, trials_clust_al] = if_align_clusters(coeffs, scores, clust_params_cell, clust_params_tr, plot_stuff);

ens_out.cells = cells_clust_al;
ens_out.trials = trials_clust_al;

%%
if plot_stuff
    f_plot_comp_scatter(coeffs(:,1:min(num_comps,3)), ens_out.cells.clust_ident);
    title(sprintf('Identified ensamble cells, %s space, %s',ensamble_method, cluster_method));
    xlabel('comp1');
    ylabel('comp2');
    zlabel('comp3');
    f_plot_comp_scatter(scores(1:min(num_comps,3),:)', ens_out.trials.clust_ident);
    title(sprintf('Identified ensamble trials, %s space, %s',ensamble_method, cluster_method));
    xlabel('comp1');
    ylabel('comp2');
    zlabel('comp3');
end


end

function clust_params = if_get_clust_params(X, clust_out)

num_dred_comps = size(X,2);

clust_label = unique(clust_out.clust_ident);
num_clust = numel(clust_label);

clust_centers = zeros(num_clust,num_dred_comps);
clust_mag = zeros(num_clust,1);
cell_num = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = clust_label(n_clust);
    clust_centers(n_clust,:) = mean(X(clust_out.clust_ident == clust1,:));
    clust_mag(n_clust) = norm(clust_centers(n_clust,:));
    cell_num(n_clust) = sum(clust_out.clust_ident == clust1);
end

[~, ens_order] = sort(clust_mag);
[~, ens_order2] = sort(ens_order);
clust_ident2 = ens_order2(clust_out.clust_ident)-1;

ens_list = cell(num_clust-1,1);
for n_ens = 1:(num_clust-1)
    ens_list{n_ens} = find(clust_ident2 == (n_ens));
end

clust_params.clust_label = clust_label - 1;
clust_params.ens_list = ens_list;
clust_params.residual_list = find(clust_ident2 == 0);
clust_params.clust_ident = clust_ident2(:);
clust_params.dend_order = clust_out.dend_order(:);
end

function [cells_clust_al, trials_clust_al] = if_align_clusters(coeffs, scores, cells_clust, trials_clust, plot_stuff)

% first align cell lists to scores sequence
ens_traces = zeros(size(scores));
for n_ens1 = 1:numel(cells_clust.ens_list)
    ens_traces(n_ens1,:) = mean(coeffs(cells_clust.ens_list{n_ens1},:)*scores);
end

best_perm_cells = f_align_comps_square(scores, ens_traces);
[~, best_perm_cells2] = sort(best_perm_cells);

ens_popl = zeros(size(coeffs));
for n_ens1 = 1:numel(trials_clust.ens_list)
    ens_popl(:,n_ens1) = mean(coeffs*scores(:,trials_clust.ens_list{n_ens1}),2);
end

best_perm_trials = f_align_comps_square(coeffs', ens_popl');
[~, best_perm_trials2] = sort(best_perm_trials);

cells_clust_al = cells_clust;
cells_clust_al.ens_list = cells_clust.ens_list(best_perm_cells);
cells_clust_al.clust_ident(cells_clust.clust_ident>0) = best_perm_cells2(cells_clust.clust_ident(cells_clust.clust_ident>0));

trials_clust_al = trials_clust;
trials_clust_al.ens_list = trials_clust.ens_list(best_perm_trials);
trials_clust_al.clust_ident(trials_clust.clust_ident>0) = best_perm_trials2(trials_clust.clust_ident(trials_clust.clust_ident>0));
end

