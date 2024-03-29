%% Hierarchical clustering demo

clear;
close all;

laod1 = load('firing_rates_cont_A1_dset1_ammn.mat');
firing_rate = laod1.firing_rate;

addpath([pwd '\functions\'])
%% cell - cell similarity

% should do clustering with normalized data, and then can plot the raster of original
est_params.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
firing_rate_norm = f_normalize(firing_rate, est_params.normalize);

hc_params = struct();
hc_params.title_tag = 'Coeffs (cells)';

hc_params.method = 'average' ; % **ward(inner square)**, **average**, single(shortest)
hc_params.distance_metric = 'cosine'; % none, euclidean, squaredeuclidean, **cosine**, hammilarity
hc_params.plot_dist_mat = 1;
hc_params.plot_clusters = 0;
hc_params.XY_label = 'cells';
hclust_out_cell = f_hcluster_wrap(firing_rate_norm, hc_params);

f_plot_raster_mean(firing_rate(hclust_out_cell.dend_order,:), 1);
title('raster cell sorted');

%% trial trial similarity
% hc_params.XY_label = 'trials';
% hclust_out_trial = f_hcluster_wrap(firing_rate', hc_params);
