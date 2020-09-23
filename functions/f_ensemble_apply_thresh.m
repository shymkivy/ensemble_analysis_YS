function ens_out = f_ensemble_apply_thresh(coeffs, scores, thresh_coeffs, thresh_scores, num_ens)

%% first detect cells

two_sided = logical(size(thresh_coeffs,2)-1);

ens_cells1 = cell(num_ens, two_sided+1);
ens_trials1 = cell(num_ens, two_sided+1);

for n_ens = 1:num_ens
    factors1 = coeffs(:,n_ens);
    [fac_sort, indx1] = sort(factors1);
    
    ens_cells1{n_ens,1} = sort(indx1(fac_sort>max(thresh_coeffs(n_ens,:))), 'ascend');
    if two_sided
        ens_cells1{n_ens,2} = sort(indx1(fac_sort<min(thresh_coeffs(n_ens,:))), 'ascend');
    end
    
    factors1 = scores(n_ens,:)';
    [fac_sort, indx1] = sort(factors1);
    ens_trials1{n_ens,1} = sort(indx1(fac_sort>max(thresh_scores(n_ens,:))), 'ascend');
    if two_sided
        ens_trials1{n_ens,2} = sort(indx1(fac_sort<min(thresh_scores(n_ens,:))), 'ascend');
    end
end
ens_cells2 = ens_cells1(:);
ens_trials2 = ens_trials1(:);

rem_ens = false(numel(ens_cells2),1);
for n_ens = 1:numel(ens_cells2)
    if numel(ens_cells2{n_ens})<2
        rem_ens(n_ens) = 1;
    end
end
ens_cells2(rem_ens) = [];
ens_trials2(rem_ens) = [];

%all_ens_cells = cat(1,ens_cells2{:});
clust_ident_cell = zeros(size(coeffs,1),1);
clust_ident_tr = zeros(size(scores,2),1);
for n_ens = 1:numel(ens_cells2)
    clust_ident_cell(ens_cells2{n_ens}) = n_ens;
    clust_ident_tr(ens_trials2{n_ens}) = n_ens;
end

% get sorting index for aethetics
LR_raster = coeffs*scores;
% sort raster so looks better
warning('off','all');
[dend_order_cell, ~] = f_hcluster(LR_raster, 'cosine', 1);
[dend_order_tr, ~] = f_hcluster(LR_raster', 'cosine', 1);
[~, inx_cell] = sort(clust_ident_cell(dend_order_cell));
[~, inx_tr] = sort(clust_ident_tr(dend_order_tr));

cells.clust_label = unique(clust_ident_cell);
cells.ens_list = ens_cells2;
cells.clust_ident = clust_ident_cell;
cells.dend_order = dend_order_cell(inx_cell);
cells.residual_list = find(clust_ident_cell == 0);

trials.clust_label = unique(clust_ident_tr);
trials.ens_list = ens_trials2;
trials.clust_ident = clust_ident_tr;
trials.dend_order = dend_order_tr(inx_tr);
trials.residual_list = find(clust_ident_tr == 0);

ens_out.trials = trials;
ens_out.cells = cells;

end

