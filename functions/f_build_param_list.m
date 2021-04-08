function params_list_full = f_build_param_list(param_struct, params_to_vary)

param_fields = fields(param_struct);
param_exists = false(numel(params_to_vary),1);
for n_param = 1:numel(params_to_vary)
    if sum(strcmpi(params_to_vary{n_param}, param_fields))
        param_exists(n_param) = 1;
    end
end

params_to_vary2 = params_to_vary(param_exists);

dims1 = ones(numel(params_to_vary2),1);
for n_param = 1:numel(params_to_vary2)
    dims1(n_param) = numel(param_struct.(params_to_vary2{n_param}));
end

num_el = prod(dims1);

params_list = cell(num_el, numel(dims1));
for n_param = 1:numel(params_to_vary2)
    if n_param == 1
        m_rep = 1;
    else
        m_rep = prod(dims1(1:(n_param-1)));
    end

    temp1 = repmat(num2cell(param_struct.(params_to_vary2{n_param})(:)'),m_rep,prod(dims1(n_param+1:end)));
    params_list(:,n_param) = temp1(:);
end

params_list_full = struct;
for n_el = 1:num_el
    for n_fl = 1:numel(param_fields)
        edit_field_ind = strcmpi(param_fields{n_fl}, params_to_vary);
        if sum(edit_field_ind)
            params_list_full(n_el).(param_fields{n_fl}) = params_list{n_el,edit_field_ind};
        else
            params_list_full(n_el).(param_fields{n_fl}) = param_struct.(param_fields{n_fl});
        end
    end
end

end