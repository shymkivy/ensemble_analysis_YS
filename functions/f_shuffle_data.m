function data_shuff = f_shuffle_data(data, shuffle_method)

if ~exist('shuffle_method', 'var')
    shuffle_method = 'circ_shift';
end
[num_cells, num_frames] = size(data);
    
data_shuff = zeros(num_cells, num_frames);
if strcmp(shuffle_method, 'circ_shift')
    for n_cell = 1:num_cells
        data_shuff(n_cell,:) = circshift(data(n_cell,:),ceil(rand(1)*num_frames));
    end
elseif strcmp(shuffle_method, 'scramble')
    % or do also scramble shuffle
    for n_cell = 1:num_cells
        data_shuff(n_cell,:) = data(n_cell, randperm(num_frames,num_frames));
    end
end

end