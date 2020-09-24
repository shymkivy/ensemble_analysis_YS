function clust_out = f_get_clust_params(X, clust_in)

num_dred_comps = size(X,2);
clust_ident = clust_in.clust_ident;
clust_label = unique(clust_ident);
num_clust = numel(clust_label);

clust_centers = zeros(num_clust,num_dred_comps);
clust_mag = zeros(num_clust,1);
num_cells = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = clust_label(n_clust);
    clust_centers(n_clust,:) = mean(X(clust_ident == clust1,:));
    clust_mag(n_clust) = norm(clust_centers(n_clust,:));
    num_cells(n_clust) = sum(clust_ident == clust1);
end

if sum(clust_label == 0)
    clust_ident2 = clust_ident + 1;
else
    clust_ident2 = clust_ident;
end

[~, ens_order] = sort(clust_mag);
[~, ens_order2] = sort(ens_order);
clust_ident3 = ens_order2(clust_ident2)-1;

ens_list = cell(num_clust-1,1);
for n_ens = 1:(num_clust-1)
    ens_list{n_ens} = find(clust_ident3 == (n_ens));
end

clust_out = clust_in;
clust_out.clust_label = unique(clust_ident3);
clust_out.ens_list = ens_list;
clust_out.non_ens_list = find(clust_ident3 == 0);
clust_out.clust_centers = clust_centers;
clust_out.clust_mag = clust_mag(ens_order);
clust_out.num_cells = num_cells(ens_order);
clust_out.clust_ident = clust_ident3(:);
end