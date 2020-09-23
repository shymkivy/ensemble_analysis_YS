function param_out = f_get_param(struct, param_name, default_val)

if isfield(struct, param_name)
    param_out = struct.(param_name);
else
    if exist('default_val', 'var') 
        param_out = default_val;
    else
        param_out = [];
    end
end

end