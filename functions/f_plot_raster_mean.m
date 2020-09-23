function f_plot_raster_mean(raster, plot_mean, trial_ind, colors1)

if ~exist('trial_ind', 'var') || isempty(trial_ind)
    plot_tr_ind = 0;
else
    plot_tr_ind = 1;
end
if ~exist('plot_mean', 'var') || isempty(plot_mean)
    plot_mean = 0;
else
    plot_mean = 1;
end

if ~exist('colors1', 'var')
    colors1 = [];
end

figure;
if ~plot_mean
    imagesc(raster);
else
    s1 = subplot(4,1,1:3);
    imagesc(raster); axis tight;
    if plot_tr_ind
        f_plot_trial_indicator3(raster, trial_ind, 1, colors1)
    end
    s2 = subplot(4,1,4);
    plot(mean(raster), 'k');
    s2.XAxis.TickValues = [];
    linkaxes([s1,s2],'x')
    subplot(s1); axis tight;
end


end