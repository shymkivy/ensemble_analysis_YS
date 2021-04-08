%% analyzing ensembles with NMF
% input: firing_rate (cells x frames)

clear;
close all;

% laod1 = load('firing_rates_cont_A1_dset1_ammn.mat');
% firing_rate = laod1.firing_rate;
% frame_rate = 30;

laod1 = load('Raster_Norm_slice.mat');
firing_rate = laod1.Raster_Norm;
frame_rate = 10;


addpath([pwd '\functions\'])
%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
% **params** are best params

estimate_params = 1;    % do estimation?
include_shuff_version = 1;
est_params.ensamble_method = 'svd';              % options: svd, nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params.smooth_SD = 0;%50:50:250;       % range of values to estimate across    % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params.num_comp = 2:4:20;               % range of values to estimate across    
est_params.shuffle_data_chunks = 0;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
est_params.reps = 1;                   % how many repeats per param 

%%
est_params.n_rep = 1:est_params.reps;
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});
if include_shuff_version
    est_params_list_s = est_params_list;
end

%% input paramseters for ensemble analysis
% NMF ensemble detection is best with thresh extraction
% for NMF best to use norm_rms(keep values positive), otherwise can also use norm_mean_std

ens_params.ensamble_method = 'nmf'; % options: svd, **nmf**, ica     % here NMF is
ens_params.num_comp = 10;
ens_params.smooth_SD = 120; % 110 is better?
ens_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensamble_extraction = 'thresh'; %  **'thresh'(only for nmf)** 'clust'(for all)
% --- for thresh detection (only nmf)
ens_params.ensamble_extraction_thresh = 'signal_z'; % 'shuff' 'signal_z' 'signal_clust_thresh'
ens_params.signal_z_thresh = 2.5;
ens_params.shuff_thresh_percent = 95;
% --- for clust detection and general sorting 
ens_params.hcluster_method = 'average';  % ward(inner square), **average**, single(shortest)     
ens_params.hcluster_distance_metric = 'cosine';  % none, euclidean, squaredeuclidean, **cosine**, hammilarity, rbf% for low component number better euclidean, otherwise use cosine
ens_params.corr_cell_thresh_percent = 95;   % to remove cells with no significant correlations
% --- other
ens_params.plot_stuff = 0;

ens_params.vol_period = 1/frame_rate*1000;

%% remove inactive cells

active_cells = sum(firing_rate,2) > 0;
firing_rate(~active_cells,:) = [];

num_cells = size(firing_rate,1);

firing_rate = firing_rate(randperm(num_cells),:);

% shuffle
firing_rate_s = f_shuffle_data(firing_rate);

firing_rate_norm = f_normalize(firing_rate, est_params.normalize);
firing_rate_norm_s = f_normalize(firing_rate_s, est_params.normalize);
%% estimate best smoothing window
if estimate_params
    est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, ens_params.vol_period);
    [~, min_ind] = min([est_params_list.test_err]);
    fprintf('From provided range, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', est_params_list(min_ind).smooth_SD, est_params.ensamble_method, est_params_list(min_ind).num_comp);

    if include_shuff_version
        fprintf('Estimating params shuff n/%d reps: ',numel(est_params_list_s));
        %dim_corr = zeros(numel(estimate_smooth_list),1);
        for n_par = 1:numel(est_params_list_s)

            params1 = est_params_list_s(n_par);
            params1.vol_period = ens_params.vol_period;
            accuracy = f_ens_estimate_corr_dim_cv(firing_rate_norm_s, params1);

            temp_fields = fields(accuracy);
            for n_fl = 1:numel(temp_fields)
                est_params_list_s(n_par).(temp_fields{n_fl}) = accuracy.(temp_fields{n_fl});
            end
            fprintf('--%d',n_par);
        end
        fprintf('\nDone\n');
        [~, min_ind] = min([est_params_list_s.test_err]);
        fprintf('From provided range, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', est_params_list_s(min_ind).smooth_SD, est_params.ensamble_method, est_params_list(min_ind).num_comp);
    end
    f_plot_cv_error_3D(est_params_list, est_params_list_s, 'smooth_SD', 'num_comp', 'test_err');
    ax1 = gca;
    ax1.Title.String = sprintf('%s L2 error from raw, (%s)',est_params.ensamble_method, ax1.Title.String);
end

%% Smooth data
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD*frame_rate/1000);

%% extract ensambles
disp('Extracting ensambles...');
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);

%% plotting ensembles
f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
title('raster cell sorted');
f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
title('raster shuffeled cell sorted');

for n_comp = 1:numel(ens_out.cells.ens_list)
    cells1 = ens_out.cells.ens_list{n_comp};
    trials1 = ens_out.trials.ens_list{n_comp};
    scores1 = ens_out.cells.ens_scores(n_comp,:);

    f_plot_ensamble_deets(firing_rate_sm, cells1, trials1, scores1);
    title([ens_params.ensamble_method ' ensamble ' num2str(n_comp)]);
end
        
%%
disp('Done');