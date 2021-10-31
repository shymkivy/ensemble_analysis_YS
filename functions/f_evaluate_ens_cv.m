function comp_test_err = f_evaluate_ens_cv(ens_out, firing_rate_norm, params)

num_comp = 1;

ens_list = ens_out.cells.ens_list;

eval_params = params;
eval_params.num_comp = num_comp;
eval_params.ensamble_method = 'pca';

comp_test_err = zeros(numel(ens_list),1);

for n_comp = 1:numel(ens_list)
    ens_raster = firing_rate_norm(ens_list{n_comp},:);
    
    accuracy = f_ens_estimate_corr_dim_cv(ens_raster, eval_params);
    comp_test_err(n_comp) = accuracy.test_err;
end

end