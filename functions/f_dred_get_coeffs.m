function [coeffs, scores, means] = f_dred_get_coeffs(dred_factors)

method = dred_factors.method;

means = dred_factors.dred_factors.means;
if strcmpi(method, 'svd')
    coeffs = dred_factors.dred_factors.coeffs;
    scores = dred_factors.dred_factors.scores;
elseif strcmpi(method, 'nmf')
    coeffs = dred_factors.dred_factors.d_W;
    scores = dred_factors.dred_factors.d_H;
elseif strcmpi(method, 'tca')
    coeffs = dred_factors.dred_factors.t_factors.U{1};
    x1 = repmat(dred_factors.dred_factors.t_factors.U{2}', 1, 1, size(dred_factors.dred_factors.t_factors.U{3},1));
    x2 = reshape(dred_factors.dred_factors.t_factors.U{3}', size(dred_factors.dred_factors.t_factors.U{3},2),1,[]);
    scores = x1.*x2.*dred_factors.dred_factors.t_factors.lambda;
elseif strcmpi(method, 'fa')
    coeffs = dred_factors.dred_factors.L;
elseif strcmpi(method, 'gpfa')
    coeffs = dred_factors.dred_factors.estParams.C;
elseif strcmpi(method, 'ica')
    coeffs = dred_factors.dred_factors.A;
    scores = dred_factors.dred_factors.icasig;
elseif strcmpi(method, 'spca')
    coeffs = dred_factors.dred_factors.SL;
end

end