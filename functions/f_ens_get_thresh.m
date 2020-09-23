function [thresh_coeffs, thresh_scores] = f_ens_get_thresh(firing_rate_ensemb, coeffs, scores, num_ens, params)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
ensamble_extraction_thresh = f_get_param(params, 'ensamble_extraction_thresh', 'shuff'); % 'signal_z' 'shuff' 'signal_clust_thresh'
plot_stuff = f_get_param(params, 'plot_stuff', 0);
signal_z_thresh = f_get_param(params, 'signal_z_thresh', 2);
shuff_thresh_percent = f_get_param(params, 'shuff_thresh_percent', 95);

shuff_rep = 50; 

two_sided = logical(sum(coeffs(:)<0));
thresh_coeffs = zeros(num_ens,two_sided+1);
thresh_scores = zeros(num_ens,two_sided+1);

thresh_coeffs_z = zeros(num_ens,two_sided+1);
thresh_scores_z = zeros(num_ens,two_sided+1);


for n_comp = 1:num_ens
    factors1 = coeffs(:,n_comp);
    center1 = median(factors1);
    z_fac = sqrt(sum((factors1-center1).^2)/(numel(factors1)-1));
    thresh_coeffs_z(n_comp,1) = center1 + signal_z_thresh*z_fac;
    if two_sided
        thresh_coeffs_z(n_comp,2) = center1 - signal_z_thresh*z_fac;
    end
end

for n_comp = 1:num_ens
    factors1 = scores(n_comp,:);
    center1 = median(factors1);
    z_fac = sqrt(sum((factors1-center1).^2)/(numel(factors1)-1));
    thresh_scores_z(n_comp) = center1 + signal_z_thresh*z_fac;
    if two_sided
        thresh_scores_z(n_comp,2) = center1 - signal_z_thresh*z_fac;
    end
end

if strcmpi(ensamble_extraction_thresh, 'shuff')
    coeffs_shuff_all = cell(shuff_rep,1);
    scores_shuff_all = cell(shuff_rep,1);
    fprintf('shuffing rep of n/%d: ', shuff_rep)
    for n_rep = 1:shuff_rep
        firing_rate_ensemb_shuff = f_shuffle_data(firing_rate_ensemb, 'circ_shift');
        train_done = 0;
        num_fail = 0;
        while ~train_done
            try
                [dred_factors_shuff, ~] = f_dred_train2(firing_rate_ensemb_shuff, num_ens, ensamble_method, 0);
                train_done = 1;
                fprintf('%d..', n_rep);
            catch
                num_fail = num_fail + 1;
                %disp('Error train, will repeat');
            end
            if num_fail > 5
                firing_rate_ensemb_shuff = f_shuffle_data(firing_rate_ensemb, 'circ_shift');
            end
        end
        [coeffs_shuff, scores_shuff] = f_dred_get_coeffs(dred_factors_shuff);
        coeffs_shuff_all{n_rep} = coeffs_shuff;
        scores_shuff_all{n_rep} = scores_shuff';
    end
    fprintf('\nDone\n');
    coeffs_shuff_all1 = cat(1,coeffs_shuff_all{:});
    scores_shuff_all1 = cat(1,scores_shuff_all{:});
    if ~two_sided
        for n_comp = 1:num_ens
            thresh_coeffs(n_comp,1) = prctile(coeffs_shuff_all1(:,n_comp), shuff_thresh_percent);
            thresh_scores(n_comp,1) = prctile(scores_shuff_all1(:,n_comp), shuff_thresh_percent); 
        end
    else
        for n_comp = 1:num_ens
            thresh_coeffs(n_comp,:) = prctile(coeffs_shuff_all1(:,n_comp), [(100-shuff_thresh_percent)/2 shuff_thresh_percent+(100-shuff_thresh_percent)/2]);
            thresh_scores(n_comp,:) = prctile(scores_shuff_all1(:,n_comp), [(100-shuff_thresh_percent)/2 shuff_thresh_percent+(100-shuff_thresh_percent)/2]);
        end
    end
elseif strcmpi(ensamble_extraction_thresh, 'signal_clust_thresh')
    for n_comp = 1:num_ens
        factors1 = coeffs(:,n_comp);
        %figure; plot(ones(numel(factors1),1),factors1, 'o')
        [~, clust_ident, ~] = f_hcluster(factors1, 'ward', two_sided+2);
        [fac_sort, indx1] = sort(factors1, 'ascend');
        clust_edge = find(abs(diff(clust_ident(indx1))));
        for ii = 1:numel(clust_edge)
            thresh_coeffs(n_comp,ii) = mean([fac_sort(clust_edge(ii)), fac_sort(clust_edge(ii)+1)]);
        end
        
        factors1 = scores(n_comp,:)';
        [~, clust_ident, ~] = f_hcluster(factors1, 'ward', two_sided+2);
        [fac_sort, indx1] = sort(factors1, 'ascend');
        clust_edge = find(abs(diff(clust_ident(indx1))));
        for ii = 1:numel(clust_edge)
            thresh_scores(n_comp,ii) = mean([fac_sort(clust_edge(ii)), fac_sort(clust_edge(ii)+1)]);
        end
    end
else
    thresh_coeffs = thresh_coeffs_z;
    thresh_scores = thresh_scores_z;
end

if plot_stuff
    max_num_plots = 10;
    width_ratio = 5; 
    for n_comp = 1:num_ens
        plot_ind = rem(n_comp-1,max_num_plots)+1;
        if plot_ind == 1
            figure;
        end
        ax1 = subplot(ceil(max_num_plots/2),2*(width_ratio+1),((plot_ind-1)*(width_ratio+1)+1):((plot_ind-1)*(width_ratio+1)+width_ratio-1)); hold on;
        stem(coeffs(:,n_comp))
        plot(ones(numel(coeffs(:,n_comp)),1)*thresh_coeffs_z(n_comp,:), '--r');
        if ~strcmpi(ensamble_extraction_thresh, 'signal_z')
            plot(ones(numel(coeffs(:,n_comp)),1)*thresh_coeffs(n_comp,:), '--g');
        end
        axis tight;
        title(sprintf('Coeffs, comp %d', n_comp));
        xlabel('Cells')
        ax2 = subplot(ceil(max_num_plots/2),2*(width_ratio+1),((plot_ind-1)*(width_ratio+1)+width_ratio)); hold on;
        [f, x] = ksdensity(coeffs(coeffs(:,n_comp)>0,n_comp),'Bandwidth',0.5);
        plot(f, x)
        ax2.YLim = ax1.YLim;
        linkaxes([ax1,ax2],'y')
    end
    
    for n_comp = 1:num_ens
        plot_ind = rem(n_comp-1,max_num_plots)+1;
        if plot_ind == 1
            figure;
        end
        subplot(ceil(max_num_plots/2),2,plot_ind); hold on;
        pl1 = plot(scores(n_comp,:));
        pl2 = plot(ones(numel(scores(n_comp,:)),1)*thresh_scores_z(n_comp,:), '--r');
        if ~strcmpi(ensamble_extraction_thresh, 'signal_z')
            pl3 = plot(ones(numel(scores(n_comp,:)),1)*thresh_scores(n_comp,:), '--g');
        end
        axis tight;
        title(sprintf('Scores, comp %d', n_comp));
        xlabel('Trials')
    end
    if ~strcmpi(ensamble_extraction_thresh, 'signal_z')
        legend([pl1 pl2(1) pl3(1)], {'data', 'z thresh', ensamble_extraction_thresh});
    else
        legend([pl1 pl2(1)], {'data', 'z thresh'});
    end
end

%figure; ecdf(coeffs_shuff_all1(:,1))

end