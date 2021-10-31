function f_plot_ensamble_deets(firing_rate_sm, cells1, trials, scores, coeffs)

[coeffs_sort, coeff_idx] = sort(coeffs, 'descend');
cells2 = cells1(coeff_idx);

raster = firing_rate_sm(cells2,:);

[num_cells, num_row] = size(raster);

trial_mark1 = ones(1,size(raster,2),3);
trial_mark2 = zeros(1,size(raster,2),3);
for n_tr = 1:numel(trials)
    trial_mark1(1,trials(n_tr),:) = [1 .6 .6];
    trial_mark2(1,trials(n_tr),:) = [1 0 0];
end

figure;
s1 = subplot(3,11,1:10); hold on;
imagesc(raster);
ylabel('Cells');
imagesc(1:num_row,num_cells+1,trial_mark2);
s1.YAxis.Direction = 'reverse';
axis tight;
s2 = subplot(3,11,12:21); hold on;
imagesc(1:size(raster,2), [0 max(mean(raster))], trial_mark1)
plot(mean(raster), 'k');
ylim([min(mean(raster)), 1.2*max(mean(raster))]);
s2.XAxis.TickValues = [];
title('mean activity');
s3 = subplot(3,11,23:32); hold on;
imagesc(1:size(raster,2), [min(scores) max(scores)], trial_mark1)
plot(scores, 'k')
title('dim reduction score')
xlabel('Frames')
ylim([1.2*min(scores), 1.2*max(scores)]);
linkaxes([s1,s2, s3],'x')
subplot(s1); axis tight;
s4 = subplot(3,11,11);
plot(coeffs_sort, 1:numel(coeffs_sort), '-ok', 'LineWidth', 2, 'MarkerSize', 6)
s4.YAxis.Direction = 'reverse';
xlabel('coeff mag');
s4.XAxis.Limits(1) = 0;
s4.XAxis.Limits(2) = s4.XAxis.Limits(2) + s4.XAxis.Limits(2)/20;
linkaxes([s1,s4],'y');
subplot(s1); axis tight;


end