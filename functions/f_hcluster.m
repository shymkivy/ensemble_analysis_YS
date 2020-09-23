function [dend_order, clust_ident, Z] = f_hcluster(data, method, num_clust)


dend_thresh = 50;
% if ~exist('fig_hand', 'var')
%     fig_hand{1} = figure;
%     fig_hand{2} = figure;
% end

if strcmpi(method, 'ward')
    Z = linkage(data,'ward');
elseif strcmpi(method, 'cosine')
    Z = linkage(pdist(data,'cosine'),'average'); % weighted
elseif strcmpi(method, 'hammilarity')
    [~, SI_hamm] = similarity_index(data,data);
    zero_diag_mat = 1 - diag(ones(size(SI_hamm,1),1));
    dist = squareform((1 - SI_hamm).*zero_diag_mat);
    Z = linkage(dist,'average'); % weighted
elseif strcmpi(method, 'rbf')
    dist = f_rbf_kernel(data, data);
    Z = linkage(dist,'ward'); % weighted
end


dend_order =  binOrd(Z);
% f1 = figure;
% [~, ~, dend_order] = dendrogram(Z, size(data,1),'ColorThreshold',dend_thresh,'Orientation','left');
% close(f1);


if exist('num_clust', 'var')
    clust_ident = cluster(Z, 'MaxClust', num_clust);
else
    clust_ident = ones(size(data,1),1);
end





end

function c = binOrd(tree)
      p = size(tree,1);
      c = zeros(1,p+1);
      k = 0;
      recBin(p);
      
      function recBin(r)
          x = tree(r,1:2)-p-1;
          if x(1)<=0
              c(k+1) = tree(r,1);
              k = k+1;
          else
              recBin(x(1));
          end
          if x(2)<=0
              c(k+1) = tree(r,2);
              k = k+1;
          else
              recBin(x(2));
          end
      end
end