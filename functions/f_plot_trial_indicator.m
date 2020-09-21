function col_width = f_plot_trial_indicator(trial_types, dend_order, num_bins, plot_y_end, sample_types_colors)
hold on;
[num_trials, ~] = size(trial_types);

trial_order = 1:num_trials;

trial_types_sort = trial_types(dend_order);
trial_order_sort = trial_order(dend_order);

gray_cmap = repmat(linspace(1,0.2,num_trials),3,1);
color_seq_tt = zeros(1,numel(trial_types),3);
color_seq_temporal = zeros(1,numel(trial_types),3);
for n_tr = 1:num_trials
    color_seq_tt(:,n_tr,:) = sample_types_colors{trial_types_sort(n_tr)}; %ops.context_types_all_colors(trial_types_sort(n_tr) == ops.context_types_all,:,:)
    color_seq_temporal(:,n_tr,:) = gray_cmap(:,trial_order_sort(n_tr));
end
color_seq_tt = repmat(color_seq_tt, num_bins ,1, 1);
color_seq_tt = permute(reshape(color_seq_tt, [],1,3),[2 1 3]);
color_seq_temporal = repmat(color_seq_temporal, num_bins,1, 1);
color_seq_temporal = permute(reshape(color_seq_temporal, [],1,3),[2 1 3]);
col_width = ceil(plot_y_end/50);

imagesc(1:num_trials*num_bins,plot_y_end+(1:col_width),repmat(color_seq_tt,col_width,1,1));
imagesc(1:num_trials*num_bins,plot_y_end+col_width+(1:col_width),repmat(color_seq_temporal,col_width,1,1));
axis tight;
end