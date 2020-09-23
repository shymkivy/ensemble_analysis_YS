function proj_out = f_ens_project_clust(coeffs, scores, clust_params_cell, clust_params_tr, raster_norm)

[num_cells, num_frames] = size(raster_norm);

fr_LR = coeffs*scores;

scores_norm = scores./vecnorm(scores')';

scores(1,:)*scores(2,:)';

proj_scores = zeros(size(scores));
for n_ens = 1:numel(clust_params_cell.ens_list)
    proj_scores(n_ens,:) = mean(coeffs(clust_params_cell.ens_list{n_ens},:)*scores);
end

trial_traces = zeros(size(scores));
for n_ens = 1:numel(clust_params_cell.ens_list)
    trial_traces(n_ens,:) = zeros(num_frames,1);
    trial_traces(n_ens,clust_params_tr.ens_list{n_ens}) = 1;
end

dist_mat = proj_scores*trial_traces';
best_perm = f_align_comps_square(dist_mat);

figure; imagesc(dist_mat)
figure; imagesc(dist_mat(:,best_perm))

for n_ens = 15
   figure; hold on
   plot(proj_scores(n_ens,:));
   plot(trial_traces(best_perm(n_ens),:));
end



proj_scores(1,:)*proj_scores(2,:)';

figure; plot(proj_scores(1,:)); axis tight
figure; plot(proj_scores(2,:)); axis tight


figure; imagesc(raster_norm(clust_params_cell.ens_list{n_ens},:))
figure; plot(mean(raster_norm(clust_params_cell.ens_list{n_ens},:))); axis tight


figure; imagesc(coeffs(clust_params_cell.ens_list{n_ens},:)*scores)

figure; plot(mean(coeffs(clust_params_cell.ens_list{n_ens},:)*scores)); axis tight


[~, ind] = max(abs(mean(coeffs(clust_params_cell.ens_list{n_ens},:)*scores_norm)*scores_norm'))

figure; plot(scores(15,:))

clust_params_cell


figure; imagesc(fr_LR(clust_params_cell.ens_list{n_ens},:))

figure; plot(mean(fr_LR(clust_params_cell.ens_list{n_ens},:)))

end