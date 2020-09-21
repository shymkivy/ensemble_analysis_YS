function cv_groups = f_make_crossval_groups(num_trials, kFold)

num_test_trials = zeros(kFold,1);

remaining_tr = num_trials;
for n_gr = 1:(kFold - 1)
    num_test_trials(n_gr) = ceil(remaining_tr/(kFold - n_gr + 1));
    remaining_tr = remaining_tr - num_test_trials(n_gr);
end
num_test_trials(kFold) = remaining_tr;

test_trial_bool = false(num_trials, kFold);
start1 = 1;
for n_gr = 1:kFold
    end1 = num_test_trials(n_gr) + start1 - 1;
    test_trial_bool(start1:end1,n_gr) = 1;
    start1 = end1 + 1;
end

cv_groups.num_test_trials = num_test_trials;
cv_groups.test_trial_bool = test_trial_bool;

end