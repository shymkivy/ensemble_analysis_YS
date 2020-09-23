function perm_out = f_align_comps_square(dist_mat_pre)

[d1, d2] = size(dist_mat_pre);


%% greedy algorithm for alignment?
dist_mat_updated = dist_mat_pre;
dist_mat_max = dist_mat_pre;
zero_mask = ones(d1,d2);
perm_col = 1:d1;
col_complete = false(d1,1);
n_itr = 0;
while ~prod(col_complete) && n_itr<10000
    n_itr = n_itr + 1;
    [~, max_ind] = max(dist_mat_max(:));
    [row,col] = ind2sub([d1 d2],max_ind);
    if row == col
        col_complete(row) = 1;
        zero_mask(row,row) = 0;
    else
        new_val = dist_mat_updated(row,col)+dist_mat_updated(col,row);
        old_val = dist_mat_updated(row,row)+dist_mat_updated(col,col);
        if new_val > old_val
            temp_coord = [perm_col(row), perm_col(col)];
            perm_col(col) = temp_coord(1);
            perm_col(row) = temp_coord(2);
            col_complete(row) = 1;
            zero_mask(row,row) = 0;
            dist_mat_updated = dist_mat_pre(:,perm_col);
        else
            zero_mask(row, col) = 0;
        end
    end
    dist_mat_max = dist_mat_pre(:,perm_col).*zero_mask;
end
% figure; imagesc(dist_mat_pre);
% title('dist pre')
% figure; imagesc(dist_mat_updated);
% title(['dist post; ' num2str(n_itr) ' iterations']);

%final_score1 = sum(diag(dist_mat_pre(:,perm_col)));
%% random perm
% dist_mat_updated = dist_mat_pre;
% col_range = 1:d1;
% perm_col = 1:d1;
% diag_sum = zeros(1000,1);
% for n_itr = 1:1000
%     col1 = randsample(col_range,1);
%     col2 = randsample(col_range(col_range ~= col1) ,1);
%     new_val = dist_mat_updated(col1,col2)+dist_mat_updated(col2,col1);
%     old_val = dist_mat_updated(col1,col1)+dist_mat_updated(col2,col2);
%     if new_val > old_val
%         temp_coord = [perm_col(col1), perm_col(col2)];
%         perm_col(col2) = temp_coord(1);
%         perm_col(col1) = temp_coord(2);
%     end
%     dist_mat_updated = dist_mat_pre(:,perm_col);
%     diag_sum(n_itr) = sum(diag(dist_mat_pre(:,perm_col)));
% end
% figure; imagesc(dist_mat_pre(:,perm_col));
% title([num2str(n_itr) ' iterations']);
% final_score2 = sum(diag(dist_mat_pre(:,perm_col)));
% figure; plot(1:1000, diag_sum);

%%
perm_out = perm_col;

end