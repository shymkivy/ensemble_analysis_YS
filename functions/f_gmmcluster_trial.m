function clust_out = f_gmmcluster_trial(X, params)
if ~exist('params', 'var')
    params = struct();
end

metric = params.metric;
[num_cells, num_comps] = size(X);

%% estimate num clust
if ~isfield(params, 'num_clust')
    
    if num_cells > 5
        k_list = 1:6;
    else
        k_list = 1:size(X,1);
    end

    E = evalclusters(X,'gmdistribution','silhouette','klist',k_list, 'Distance', 'sqEuclidean');
    num_clust = E.OptimalK;
    
else
    num_clust = params.num_clust;
end

%%



threshold = sqrt(chi2inv(0.99,2));

k = num_clust; % Number of GMM components
options = statset('MaxIter',1000);
cov_type = 'full';
gmfit = fitgmdist(X,k,'RegularizationValue',params.RegularizationValue,'CovarianceType',cov_type, ...
    'SharedCovariance',false,'Options',options); % Fitted GMM
clusterX = cluster(gmfit,X); % Cluster index 
if num_comps == 2
    d = 500; % Grid length
    x1 = linspace(min(X(:,1))-2, max(X(:,2))+2, d);
    x2 = linspace(min(X(:,1))-2, max(X(:,2))+2, d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];


    mahalDist = mahal(gmfit,X0); % Distance from each grid point to each GMM component
    % Draw ellipsoids over each GMM component and show clustering result.
    figure;
    h1 = gscatter(X(:,1),X(:,2),clusterX);
    hold on
        for m = 1:k
            idx = mahalDist(:,m)<=threshold;
            Color = h1(m).Color*0.75 - 0.5*(h1(m).Color - 1);
            h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
            uistack(h2,'bottom');
        end    
    plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
    title(sprintf('Sigma is %s\nSharedCovariance = %s',cov_type,'false'),'FontSize',8)
    legend(h1,{'1','2','3'})
    hold off;
end

clust_out.num_clust = num_clust;
clust_out.clust_ident = clusterX;
clust_out.gmfit = gmfit;


end