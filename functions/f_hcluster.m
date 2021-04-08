function [dend_order, clust_ident, Z] = f_hcluster(data, method, num_clust)

Z = linkage(data,method);

dend_order =  binOrd(Z)';

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