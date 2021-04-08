function est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, vol_period)

fprintf('Estimating params n/%d reps: ',numel(est_params_list));
%dim_corr = zeros(numel(estimate_smooth_list),1);
for n_par = 1:numel(est_params_list)

    params1 = est_params_list(n_par);
    params1.vol_period = vol_period;
    accuracy = f_ens_estimate_corr_dim_cv(firing_rate_norm, params1);

    temp_fields = fields(accuracy);
    for n_fl = 1:numel(temp_fields)
        est_params_list(n_par).(temp_fields{n_fl}) = accuracy.(temp_fields{n_fl});
    end
    fprintf('--%d',n_par);
end
fprintf('\nDone\n');
[~, min_ind] = min([est_params_list.test_err]);

end