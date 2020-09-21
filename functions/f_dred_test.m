function Ycs = f_dred_test(data_3d, dred_factors, method)

[num_cells, num_bins, num_trials] = size(data_3d);
data_2d = reshape(data_3d, num_cells, []);

%%
Ycs = zeros(size(data_2d)); 
d = dred_factors.means;
if strcmpi(method, 'svd')
    L = dred_factors.coeffs;
    for n_cell = 1:num_cells  
        % Indices 1:yDim with i removed
        mi = [1:(n_cell-1) (n_cell+1):num_cells];
        % inv(A)*B ~ A\B  and B*inv(A) ~ B/A
        % inv(L'*L) * L'*(Y-mean)
        Xmi = (L(mi,:)' * L(mi,:)) \ L(mi,:)' * bsxfun(@minus, data_2d(mi,:), d(mi));
        Ycs(n_cell,:) = L(n_cell,:) * Xmi + d(n_cell);  
    end
elseif strcmpi(method, 'nmf')
    d_W = dred_factors.d_W;
    for n_cell = 1:num_cells  
        % Indices 1:yDim with i removed
        mi = [1:(n_cell-1) (n_cell+1):num_cells];
        % inv(A)*B ~ A\B  and B*inv(A) ~ B/A
        % inv(L'*L) * L'*(Y-mean)
        Xmi = (d_W(mi,:)' * d_W(mi,:)) \ d_W(mi,:)' * data_2d(mi,:);
        Ycs(n_cell,:) = d_W(n_cell,:) * Xmi;  
    end
elseif strcmpi(method, 'fa')
    [Ycs, ~] = cosmoother_fa(data_2d, dred_factors);
elseif strcmpi(method, 'tca')
    
    Ycs_3d = reshape(Ycs,num_cells,num_bins,[]);
    tf = dred_factors.t_factors;
    
    tca_met = 2;
    if tca_met == 1
        % this fits both cells and traces and computes trial variables
        num_rank = numel(tf.lambda);   
        L = zeros(num_rank,num_cells,size(tf.U{2},1));
        for n_rank = 1:num_rank
            L(n_rank,:,:) = tf.lambda(n_rank)*tf.U{1}(:,n_rank)* tf.U{2}(:,n_rank)';
        end
        for n_cell = 1:num_cells  
            % Indices 1:yDim with i removed
            mi = [1:(n_cell-1) (n_cell+1):num_cells];
            LLmi = double(full(ttt(tensor(L(:,mi,:)),tensor(L(:,mi,:)),[2 3])));
            LY = double(full(ttt(tensor(L(:,mi,:)), tensor(data_3d(mi,:,:)-d(mi)), [2 3], [1 2])));
            Xmi = LLmi \ LY;
            Ycs_3d(n_cell,:,:) = ttt(tensor(L(:,n_cell,:)),tensor(Xmi),1,1)+d(n_cell);  
            % inv(A)*B ~ A\B  and B*inv(A) ~ B/A
            % inv(L'*L) * L'*(Y-mean)
            %Xmi = (L(mi,:)' * L(mi,:)) \ L(mi,:)' * data_2d(mi,:);
            %Ycs(n_cell,:) = d_W(n_cell,:) * Xmi;  
        end
    elseif tca_met == 2
        % this only fits cells and guesses traces + trial vars (more stringent)
        L = tf.lambda'.* tf.U{1};
        for n_cell = 1:num_cells  
            % Indices 1:yDim with i removed
            mi = [1:(n_cell-1) (n_cell+1):num_cells];
            iLLmi = inv(L(mi,:)'*L(mi,:));
            LY = ttt(tensor(L(mi,:)'), tensor(data_3d(mi,:,:)-d(mi)), 2, 1);
            Xmi = ttt(tensor(iLLmi),LY,2,1);
            Ycs_3d(n_cell,:,:) = ttt(tensor(L(n_cell,:)), Xmi,2,1)+d(n_cell);  
            % inv(A)*B ~ A\B  and B*inv(A) ~ B/A
            % inv(L'*L) * L'*(Y-mean)
            %Xmi = (L(mi,:)' * L(mi,:)) \ L(mi,:)' * data_2d(mi,:);
            %Ycs(n_cell,:) = d_W(n_cell,:) * Xmi;  
        end
    end

    Ycs = reshape(Ycs_3d, num_cells,[]);
    %data__sn3d = reshape(data_n2d, num_cells, num_bins, num_trials);
    %est_factors = cp_als(tensor(data__sn3d),num_comp);
    %dred_factors = est_factors;
    %dred_data = full(est_factors) + data_means;
    %cp_wopt can optimize with data held out
    
elseif strcmpi(method, 'gpfa')
    seq = struct();
    for n_tr = 1:num_trials
        seq(n_tr).trialId = n_tr;
        seq(n_tr).T = num_bins;
        seq(n_tr).segId = 1;
        seq(n_tr).y = squeeze(data_3d(:,:,n_tr));
    end
    estParams = dred_factors.estParams;
    num_comp = size(estParams.C,2);
    % Leave-neuron-out prediction on test data 
    if estParams.notes.RforceDiagonal
      seq = cosmoother_gpfa_viaOrth_fast(seq, estParams, 1:num_comp);
    else
      seq = cosmoother_gpfa_viaOrth(seq, estParams, 1:num_comp);
    end
    
    fn = sprintf('ycsOrth%.2d',num_comp);
    Ycs = ([seq.(fn)]);
end



end