function data_smooth = f_smooth_gauss(data, smoothSDbinSizeRatio)

% smoothSDbinRatio is the sigma in bins

[num_cells, num_bins, num_trials] = size(data);

%% Smooth data
% compute the size of smoothing window

if smoothSDbinSizeRatio>0
    kernel_half_size = ceil(sqrt(-log(0.05)*2*smoothSDbinSizeRatio^2));
    gaus_win = -kernel_half_size:kernel_half_size;
    gaus_kernel = exp(-((gaus_win).^2)/(2*smoothSDbinSizeRatio^2));
    gaus_kernel = gaus_kernel/sum(gaus_kernel);

    %figure; plot(gaus_win,gaus_kernel)

    data_smooth = zeros(num_cells, num_bins, num_trials);
    for n_cell = 1:num_cells
        for n_trial = 1:num_trials
            data_smooth(n_cell,:,n_trial) = conv(data(n_cell,:,n_trial), gaus_kernel, 'same');
        end
    end
else
    data_smooth = data;
end

end