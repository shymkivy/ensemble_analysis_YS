function D = f_pdist2_YS(X, Y, distance_metric)
if ~exist('distance_metric', 'var') || isempty(distance_metric)
    distance_metric = 'euclidean';
end

if strcmpi(distance_metric, 'hammilarity')
    [~, SI_hamm] = similarity_index(X,Y);
    D = (1 - SI_hamm);
elseif strcmpi(distance_metric, 'rbf')
    D = 1 - f_rbf_kernel(X, Y);
else
    D = pdist2(X, Y, distance_metric);
end

end