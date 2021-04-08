function f_plot_ensamble_deets(firing_rate_sm, cells1, trials, scores1)

raster = firing_rate_sm(cells1,:);

[num_cells, num_row] = size(raster);

trial_mark1 = ones(1,size(raster,2),3);
trial_mark2 = zeros(1,size(raster,2),3);
for n_tr = 1:numel(trials)
    trial_mark1(1,trials(n_tr),:) = [1 .6 .6];
    trial_mark2(1,trials(n_tr),:) = [1 0 0];
end

figure;
s1 = subplot(3,1,1); hold on;
imagesc(raster);
ylabel('Cells');
imagesc(1:num_row,num_cells+1,trial_mark2);
s1.YAxis.Direction = 'reverse';
axis tight;
s2 = subplot(3,1,2); hold on;
imagesc(1:size(raster,2), [0 max(mean(raster))], trial_mark1)
plot(mean(raster), 'k');
ylim([min(mean(raster)), 1.2*max(mean(raster))]);
s2.XAxis.TickValues = [];
title('mean activity');
s3 = subplot(3,1,3); hold on;
imagesc(1:size(raster,2), [min(scores1) max(scores1)], trial_mark1)
plot(scores1, 'k')
title('dim reduction score')
xlabel('Frames')
ylim([1.2*min(scores1), 1.2*max(scores1)]);
linkaxes([s1,s2, s3],'x')
subplot(s1); axis tight;


end