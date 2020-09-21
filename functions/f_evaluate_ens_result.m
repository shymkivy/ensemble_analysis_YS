function data_out = f_evaluate_ens_result(ens_output, ground_truth, plot_stuff)

[gt_mat, gt_list] = if_list_2_mat(ground_truth, []);
[data_ens_mat, ens_list] = if_list_2_mat(ens_output, []);


SI_gt_data = 1 - pdist2(gt_mat,data_ens_mat,'cosine');%similarity_index(gt_mat, data_ens_mat);

if plot_stuff
    figure;
    subplot(2,2,1);
    imagesc(SI_gt_data)
    ylabel('ground truth ens');
    xlabel('data ens');
    title('clust pre');
end
%[~, prem_ind] = max(SI_gt_data);

%SI_gt_data2 = SI_gt_data;
%sizSI = size(SI_gt_data2);
% clust_prem_ind = zeros(numel(gt_list),1);
% for n_tt = 1:(numel(gt_list)-1)
%     [~, temp_ind] = max(SI_gt_data2(:));
%     [gt_ind, dt_ind] = ind2sub(sizSI, temp_ind);
%     clust_prem_ind(gt_ind) = dt_ind;
%     SI_gt_data2(:,dt_ind) = 0;
% end
%clust_prem_ind(clust_prem_ind==0) = find(~sum(ens_list' == clust_prem_ind(clust_prem_ind>0)));

max_size = max([numel(gt_list), numel(ens_list)]);
all_perms = perms(1:max_size);
num_perms = size(all_perms,1);
acc1 = zeros(num_perms,1);
for n_pr = 1:size(all_perms,1)
    if numel(ens_list)>= numel(gt_list)
        temp_SI = SI_gt_data(:,all_perms(n_pr,:));
        perm_ens = 1;
    else
        temp_SI = SI_gt_data(all_perms(n_pr,:),:);
        perm_ens = 0;
    end
    if min(size(temp_SI)) == 1
        acc1(n_pr) = temp_SI(1,1);
    else
        acc1(n_pr) = sum(diag(temp_SI));
    end
    
end
[~, max_perm_ind] = max(acc1);

best_perm = all_perms(max_perm_ind,:);

if perm_ens
    data_ens_mat_sort = data_ens_mat(best_perm,:);
    gt_mat_sort = gt_mat;
else
    data_ens_mat_sort = data_ens_mat;
    gt_mat_sort = gt_mat(best_perm,:);
end

SI_gt_data3 = 1 - pdist2(gt_mat_sort,data_ens_mat_sort,'cosine');

core_size = min([numel(gt_list), numel(ens_list)]);
if numel(ens_list)~= numel(gt_list)
    if perm_ens
        [acc_vals, ens_num1] = max(SI_gt_data3(:,core_size+1:end),[],1);
    else
        [acc_vals, ens_num1] = max(SI_gt_data3(core_size+1:end,:),[],2);
    end
else
    ens_num1 = [];
    acc_vals = [];
end

if perm_ens
   gt_seq = [gt_list; ens_num1];
   ens_seq = ens_list(best_perm);
else
   gt_seq = gt_list(best_perm);
   ens_seq = [ens_list; ens_num1];
end

[~, ind1] = sort(ens_seq);
aligned_seq = [gt_seq, ens_seq];
aligned_seq = aligned_seq(ind1,:);


if plot_stuff
    subplot(2,2,2);
    imagesc(SI_gt_data3);
    ylabel('ground truth ens');
    xlabel('data ens');
    title('clust aligned');
    
    if core_size > 1
        bar_plot1 = diag(SI_gt_data3(1:core_size, 1:core_size))';
    else
        bar_plot1 = SI_gt_data3(1:core_size, 1:core_size);
    end
    subplot(2,2,3);
    bar([bar_plot1'; acc_vals(:)]); axis tight;
    ylim([0 1]);
    xlabel('cluster num');
    ylabel('overlap accuracy');
    title('detection accuracy');
    
    num_gt = zeros(max_size,1);
    for ii = 1:numel(gt_list)
        num_gt(ii) = sum(gt_mat(ii,:),2);
    end
    num_ens = zeros(max_size,1);
    for ii = 1:numel(ens_list)
        num_ens(ii) = sum(data_ens_mat_sort(ii,:),2);
    end
    
    subplot(2,2,4);
    bar([num_gt, num_ens]); axis tight;
    xlabel('cluster num');
    ylabel('num cells');
    title('ens size');
end

if core_size > 1
    accuracy1 = diag(SI_gt_data3);
else
    accuracy1 = SI_gt_data3(1);
end
data_out.ens_perm_ind = best_perm;
data_out.accuracy = accuracy1;
data_out.ground_truth_cat_order = gt_list;
data_out.data_cat_order = ens_list;
data_out.aligned_seq = aligned_seq;

end

function [ens_mat, ens_types] = if_list_2_mat(ens_list, ens_types)

num_cells = max(cat(1,ens_list{:}));

if exist('ens_types', 'var') || isempty(ens_types)
    ens_types = (1:numel(ens_list))'-1;
end

ens_mat = zeros(numel(ens_list),num_cells);
for n_cl=1:numel(ens_types)
    ens_mat(n_cl,ens_list{n_cl}) = 1;
end

end