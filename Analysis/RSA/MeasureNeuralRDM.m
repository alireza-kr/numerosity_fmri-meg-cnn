%__________________________________________________________________________
function rdm = MeasureNeuralRDM(ds,mode,params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-NONCROSS-VALIDATED PEARSON CORRELATION WITHOUT MNN
%==========================================================================
if strcmp(mode,'Pearson-none-none-dnn')

%-Center the data
samples = ds;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end

%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'correlation'));
    
%-NONCROSS-VALIDATED SPEARMAN CORRELATION WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Spearman-none-none-dnn')

%-Center the data
samples = ds;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end
    
%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'spearman'));    

%-NONCROSS-VALIDATED EUCLIDEAN WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Euclidean-none-none-dnn')

%-Center the data
samples = ds;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end
    
%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'euclidean'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-NONCROSS-VALIDATED PEARSON CORRELATION WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Pearson-none-none-mri')

%-Center the data
samples = ds.samples;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end

%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'correlation'));
    
%-NONCROSS-VALIDATED SPEARMAN CORRELATION WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Spearman-none-none-mri')

%-Center the data
samples = ds.samples;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end
    
%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'spearman'));    

%-NONCROSS-VALIDATED EUCLIDEAN WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Euclidean-none-none-mri')

%-Center the data
samples = ds.samples;
if params.center_data
    samples = bsxfun(@minus,samples,mean(samples,1));
end
    
%-Compute the pair-wise distance between all dataset samples using cosmo_pdist
rdm = cosmo_squareform(cosmo_pdist(samples,'euclidean'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-NONCROSS-VALIDATED PEARSON CORRELATION WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'Pearson-none-none-meg')

X = ds.trial;
y = single(ds.index);

n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_conditions, n_conditions, n_time);

for c = 1:n_conditions
    sample(c,:,:) = mean(X(y==c,:,:),1);
end

for t = 1:n_time
    result(:,:,t) = cosmo_squareform(cosmo_pdist(sample(:,:,t),'correlation'));
end

rdm = result;

%-NONCROSS-VALIDATED PEARSON CORRELATION WITH MNN
%==========================================================================
elseif strcmp(mode,'Pearson-mnn-none-meg')

X = ds.trial;
y = single(ds.index);

n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

X = whiten_data(X,y);

result = nan(n_conditions, n_conditions, n_time);

for c = 1:n_conditions
    sample(c,:,:) = mean(X(y==c,:,:),1);
end

for t = 1:n_time
    result(:,:,t) = cosmo_squareform(cosmo_pdist(sample(:,:,t),'correlation'));
end

rdm = result;

%-CROSS-VALIDATED PEARSON CORRELATION WITH MNN
%==========================================================================
elseif strcmp(mode,'Pearson-mnn-cv-meg')

X = ds.trial;
y = single(ds.index);

n_perm = params.perm;      % number of permutations
n_pseudo = params.pseudo;  % number of pseudo-trials
n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_perm, n_conditions, n_conditions, n_time);

conditions = unique(y);
n_trials = histc(y, conditions);

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)

    %-precompute permutations
    ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    ind_pseudo_test = nan(n_conditions, n_conditions, 2);
    labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    labels_pseudo_test = nan(n_conditions, n_conditions, 2);
    for c1 = 1:n_conditions
        range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for c2 = 1:n_conditions
            range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
            ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
            ind_pseudo_test(c1, c2, :) = [c1 c2];
            labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                [conditions(c1)*ones(1, n_pseudo - 1) conditions(c2)*ones(1, n_pseudo - 1)];
            labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
        end
    end              
    train_indices = cell(1, n_conditions*(n_pseudo-1));
    test_indices = cell(1, n_conditions);
    %-separate permutation for each class
    for c1 = 1:n_conditions
        prm_ = randperm(n_trials(c1));                
        prm = cell(1, n_pseudo);
        splitsize = n_trials(c1) / n_pseudo;
        for i = 1:n_pseudo
            idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
            prm{i} = prm_(idxs + 1);
        end                                
        ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
        xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for i = 1:length(xrange)
            train_indices{xrange(i)} = ind{i};
        end
        test_indices{c1} = ind{end};
    end                                

    %-1. Compute pseudo-trials for training and test
    Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
    Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
    for i = 1:length(train_indices)
        Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
    end
    for i = 1:length(test_indices)
        Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
    end

    %-2. Whitening using the Epoch method
    sigma_conditions = reshape(squeeze(labels_pseudo_train(1, :, n_pseudo:end))', 1, []);
    sigma_ = nan(n_conditions, n_sensors, n_sensors);
    for c = 1:n_conditions
        %-compute sigma for each time point, then average across time
        tmp_ = nan(n_time, n_sensors, n_sensors);
        for t = 1:n_time
            tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
        end
        sigma_(c, :, :) = mean(tmp_, 1);
    end
    %-average across conditions
    sigma = squeeze(mean(sigma_, 1));
    sigma_inv = sigma^-0.5;
    for t = 1:n_time
        Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
        Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
    end

    for t = 1:n_time
        for c1 = 1:n_conditions-1
            for c2 = c1+1:n_conditions
                %-3. Apply distance measure to training data
                data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_train(c1, c2, :));
                classes = unique(y_train);
                %-Pearson
                A1_ps = mean(data_train(y_train==classes(1), :), 1);
                B1_ps = mean(data_train(y_train==classes(2), :), 1);
                var_A1_ps = var(A1_ps);
                var_B1_ps = var(B1_ps);
                denom_noncv_ps = sqrt(var_A1_ps * var_B1_ps);

                %-4. Validate distance measure on testing data
                data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                y_test = squeeze(labels_pseudo_test(c1, c2, :));

                %-Pearson
                A2_ps = mean(data_test(y_test==classes(1), :), 1);
                B2_ps = mean(data_test(y_test==classes(2), :), 1);
                cov_a1b2_ps = getfield(cov(A1_ps, B2_ps), {2});
                cov_b1a2_ps = getfield(cov(B1_ps, A2_ps), {2});
                cov_ab_ps = (cov_a1b2_ps + cov_b1a2_ps) / 2;
                var_A12_ps = getfield(cov(A1_ps, A2_ps), {2});
                var_B12_ps = getfield(cov(B1_ps, B2_ps), {2});
                reg_factor_var = 0.1; reg_factor_denom = 0.25;   % regularization
                denom_ps = sqrt(max(reg_factor_var * var_A1_ps, var_A12_ps) * max(reg_factor_var * var_B1_ps, var_B12_ps));
                denom_ps = max(reg_factor_denom * denom_noncv_ps, denom_ps);
                r_ps = cov_ab_ps / denom_ps; 
                r_ps = min(max(-1, r_ps), 1);
                result(f, c1, c2, t) = 1 - r_ps;
                result(f, c2, c1, t) = result(f, c1, c2, t);
            end
        end
    end
end

%-average across permutations
result = squeeze(nanmean(result, 1));

for c = 1:n_conditions
    result(c,c,:) = 0;
end

rdm = result;

%-NONCROSS-VALIDATED EUCLIDEAN WITH MNN
%==========================================================================
%-Squared Euclidean distance on multivariately noise-normalized patterns
%is equivalent to the squared Mahalanobis distance
%https://www.sciencedirect.com/science/article/pii/S1053811915011258
%https://www.sciencedirect.com/science/article/pii/S1053811918301411
elseif strcmp(mode,'Euclidean-mnn-none-meg')

X = ds.trial;
y = single(ds.index);

n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

X = whiten_data(X,y);

result = nan(n_conditions, n_conditions, n_time);

for c = 1:n_conditions
    sample(c,:,:) = mean(X(y==c,:,:),1);
end

for t = 1:n_time
    result(:,:,t) = cosmo_squareform(cosmo_pdist(sample(:,:,t),'euclidean'));
end

rdm = result;

%-CROSS-VALIDATED EUCLIDEAN WITH MNN
%==========================================================================
% Cross-validated mahalanobis distance is equivalent to
%linear discriminant contrast (LDC)
%https://www.sciencedirect.com/science/article/pii/S1053811915011258
%https://www.sciencedirect.com/science/article/pii/S1053811918301411
elseif strcmp(mode,'Euclidean-mnn-cv-meg')

X = ds.trial;
y = single(ds.index);

n_perm = params.perm;      % number of permutations
n_pseudo = params.pseudo;  % number of pseudo-trials
n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_perm, n_conditions, n_conditions, n_time);

conditions = unique(y);
n_trials = histc(y, conditions);

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)

    %-precompute permutations
    ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    ind_pseudo_test = nan(n_conditions, n_conditions, 2);
    labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    labels_pseudo_test = nan(n_conditions, n_conditions, 2);
    for c1 = 1:n_conditions
        range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for c2 = 1:n_conditions
            range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
            ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
            ind_pseudo_test(c1, c2, :) = [c1 c2];
            labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                [conditions(c1)*ones(1, n_pseudo - 1) conditions(c2)*ones(1, n_pseudo - 1)];
            labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
        end
    end              
    train_indices = cell(1, n_conditions*(n_pseudo-1));
    test_indices = cell(1, n_conditions);
    %-separate permutation for each class
    for c1 = 1:n_conditions
        prm_ = randperm(n_trials(c1));                
        prm = cell(1, n_pseudo);
        splitsize = n_trials(c1) / n_pseudo;
        for i = 1:n_pseudo
            idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
            prm{i} = prm_(idxs + 1);
        end                                
        ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
        xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for i = 1:length(xrange)
            train_indices{xrange(i)} = ind{i};
        end
        test_indices{c1} = ind{end};
    end                                

    %-1. Compute pseudo-trials for training and test
    Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
    Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
    for i = 1:length(train_indices)
        Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
    end
    for i = 1:length(test_indices)
        Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
    end

    %-2. Whitening using the Epoch method
    sigma_conditions = reshape(squeeze(labels_pseudo_train(1, :, n_pseudo:end))', 1, []);
    sigma_ = nan(n_conditions, n_sensors, n_sensors);
    for c = 1:n_conditions
        %-compute sigma for each time point, then average across time
        tmp_ = nan(n_time, n_sensors, n_sensors);
        for t = 1:n_time
            tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
        end
        sigma_(c, :, :) = mean(tmp_, 1);
    end
    %-average across conditions
    sigma = squeeze(mean(sigma_, 1));
    sigma_inv = sigma^-0.5;
    for t = 1:n_time
        Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
        Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
    end

    for t = 1:n_time
        for c1 = 1:n_conditions-1
            for c2 = c1+1:n_conditions
                %-3. Apply distance measure to training data
                data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_train(c1, c2, :));
                classes = unique(y_train);
                %-Euclidean
                dist_train_ec = mean(data_train(y_train==classes(1), :), 1) - ...
                    mean(data_train(y_train==classes(2), :), 1);

                %-4. Validate distance measure on testing data
                data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                y_test = squeeze(labels_pseudo_test(c1, c2, :));

                %-Euclidean
                dist_test_ec = mean(data_test(y_test==classes(1), :), 1) - ...
                    mean(data_test(y_test==classes(2), :), 1);
                result(f, c1, c2, t) = dot(dist_train_ec, dist_test_ec);
                result(f, c2, c1, t) = result(f, c1, c2, t);
            end
        end
    end
end

%-average across permutations
result = squeeze(nanmean(result, 1));

for c = 1:n_conditions
    result(c,c,:) = 0;
end

rdm = result;

%-SVM (v1) WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'SVM.v1-none-none-meg')

X = ds.trial;
y = single(ds.index);

n_perm = params.perm;      % number of permutations
n_pseudo = params.pseudo;  % number of pseudo-trials
n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_perm, n_conditions, n_conditions, n_time);

conditions = unique(y);
n_trials = histc(y, conditions);

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)

    %-precompute permutations
    ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    ind_pseudo_test = nan(n_conditions, n_conditions, 2);
    labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    labels_pseudo_test = nan(n_conditions, n_conditions, 2);
    for c1 = 1:n_conditions
        range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for c2 = 1:n_conditions
            range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
            ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
            ind_pseudo_test(c1, c2, :) = [c1 c2];
            labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                [conditions(c1)*ones(1, n_pseudo - 1) conditions(c2)*ones(1, n_pseudo - 1)];
            labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
        end
    end              
    train_indices = cell(1, n_conditions*(n_pseudo-1));
    test_indices = cell(1, n_conditions);
    %-separate permutation for each class
    for c1 = 1:n_conditions
        prm_ = randperm(n_trials(c1));                
        prm = cell(1, n_pseudo);
        splitsize = n_trials(c1) / n_pseudo;
        for i = 1:n_pseudo
            idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
            prm{i} = prm_(idxs + 1);
        end                                
        ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
        xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for i = 1:length(xrange)
            train_indices{xrange(i)} = ind{i};
        end
        test_indices{c1} = ind{end};
    end                                

    %-1. Compute pseudo-trials for training and test
    Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
    Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
    for i = 1:length(train_indices)
        Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
    end
    for i = 1:length(test_indices)
        Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
    end

    for t = 1:n_time
        for c1 = 1:n_conditions-1
            for c2 = c1+1:n_conditions                        
                %-2. Fit the classifier using training data
                data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_train(c1, c2, :));
                model_svm = svmtrain(y_train, data_train, '-t 0 -q 0');

                %-3. Compute predictions using test data
                data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_test(c1, c2, :));
                predictions = svmpredict(y_train, data_test, model_svm, '-q 0');

                %-4. Compute dissimilarity and store in RDM
                dissimilarity = mean(predictions == y_train) - 0.5;
                result(f, c1, c2, t) = dissimilarity;
                result(f, c2, c1, t) = result(f, c1, c2, t);
            end
        end
    end
end

%-average across permutations
result = squeeze(nanmean(result, 1));

for c = 1:n_conditions
    result(c,c,:) = 0;
end

rdm = result;

%-SVM (v1) WITH MNN
%==========================================================================
elseif strcmp(mode,'SVM.v1-mnn-none-meg')

X = ds.trial;
y = single(ds.index);

n_perm = params.perm;      % number of permutations
n_pseudo = params.pseudo;  % number of pseudo-trials
n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_perm, n_conditions, n_conditions, n_time);

conditions = unique(y);
n_trials = histc(y, conditions);

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)

    %-precompute permutations
    ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    ind_pseudo_test = nan(n_conditions, n_conditions, 2);
    labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    labels_pseudo_test = nan(n_conditions, n_conditions, 2);
    for c1 = 1:n_conditions
        range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for c2 = 1:n_conditions
            range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
            ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
            ind_pseudo_test(c1, c2, :) = [c1 c2];
            labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                [conditions(c1)*ones(1, n_pseudo - 1) conditions(c2)*ones(1, n_pseudo - 1)];
            labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
        end
    end              
    train_indices = cell(1, n_conditions*(n_pseudo-1));
    test_indices = cell(1, n_conditions);
    %-separate permutation for each class
    for c1 = 1:n_conditions
        prm_ = randperm(n_trials(c1));                
        prm = cell(1, n_pseudo);
        splitsize = n_trials(c1) / n_pseudo;
        for i = 1:n_pseudo
            idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
            prm{i} = prm_(idxs + 1);
        end                                
        ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
        xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for i = 1:length(xrange)
            train_indices{xrange(i)} = ind{i};
        end
        test_indices{c1} = ind{end};
    end                                

    %-1. Compute pseudo-trials for training and test
    Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
    Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
    for i = 1:length(train_indices)
        Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
    end
    for i = 1:length(test_indices)
        Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
    end

    %-2. Whitening using the Epoch method
    sigma_conditions = reshape(squeeze(labels_pseudo_train(1, :, n_pseudo:end))', 1, []);
    sigma_ = nan(n_conditions, n_sensors, n_sensors);
    for c = 1:n_conditions
        %-compute sigma for each time point, then average across time
        tmp_ = nan(n_time, n_sensors, n_sensors);
        for t = 1:n_time
            tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
        end
        sigma_(c, :, :) = mean(tmp_, 1);
    end
    %-average across conditions
    sigma = squeeze(mean(sigma_, 1));
    sigma_inv = sigma^-0.5;
    for t = 1:n_time
        Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
        Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
    end

    for t = 1:n_time
        for c1 = 1:n_conditions-1
            for c2 = c1+1:n_conditions                        
                %-3. Fit the classifier using training data
                data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_train(c1, c2, :));
                model_svm = svmtrain(y_train, data_train, '-c 1 -q 0 -t 0');

                %-4. Compute predictions using test data
                data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_test(c1, c2, :));
                predictions = svmpredict(y_train, data_test, model_svm, '-q 0 -t 0');

                %-5. Compute dissimilarity and store in RDM
                dissimilarity = mean(predictions == y_train) - 0.5;
                result(f, c1, c2, t) = dissimilarity;
                result(f, c2, c1, t) = result(f, c1, c2, t);
            end
        end
    end
end

%-average across permutations
result = squeeze(nanmean(result, 1));

for c = 1:n_conditions
    result(c,c,:) = 0;
end

rdm = result;

%-SVM (v2) WITHOUT MNN
%==========================================================================
elseif strcmp(mode,'SVM.v2-none-cv-meg')

%-@Martin Hebart
%-X: trial x channels x time
%-Column vector called label of length pseudotrial with conditions 1 to 32
%-Column vector called chunk (of length pseudotrial with index within condition (e.g. 1 to 5))

X = ds.trial;
y = single(ds.index);

n_perm = params.perm;      % number of permutations
n_pseudo = params.pseudo;  % number of pseudo-trials
n_conditions = length(unique(y));
n_trials = size(X, 1);
n_sensors = size(X, 2);
n_time = size(X, 3);

result = nan(n_perm, n_conditions, n_conditions, n_time);

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)

    %-Create pseudo trials
    %..........................................................................
    measure_args = struct();
    measure_args.ratio = 1/n_pseudo;
    measure_args.resamplings = 1;
    measure_args.split_by = {'targets'};
    %https://it.mathworks.com/matlabcentral/answers/116502-concatenation-of-3d-array-into-2d-array
    ds_pseudo.samples = X(:,:);
    ds_pseudo.sa.targets = y';
    ds_pseudo = cosmo_average_samples(ds_pseudo,measure_args);
    label = ds_pseudo.sa.targets;
    %..........................................................................

    %-Chunkize the data
    %..........................................................................
    ds_chunk.samples = ds_pseudo.samples;
    ds_chunk.sa.targets = ds_pseudo.sa.targets;
    ds_chunk.sa.chunks=(1:size(ds_pseudo.samples,1))';
    chunk = cosmo_chunkize(ds_chunk,nnz(ds_pseudo.sa.targets==1));
    %..........................................................................

    %-Make the design
    %..........................................................................
    cfg.files.label = double(label);
    cfg.files.chunk = double(chunk);
    design = make_design_cv(cfg);
    clear cfg
    %..........................................................................

    n_steps = size(design.train,2);
    n_trials = size(ds_pseudo.samples,1);
    Xpseudo = reshape(ds_pseudo.samples,[n_trials n_sensors n_time]);

    clear result;
    for i_time = n_time:-1:1 % trick to run without preallocation
        % precompute kernel to save a lot of time
        %https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
        %https://stackoverflow.com/questions/2474460/precomputed-kernels-with-libsvm-in-python
        %https://stackoverflow.com/questions/7715138/using-precomputed-kernels-with-libsvm
        kernel = Xpseudo(:,:,i_time)*Xpseudo(:,:,i_time)';

        % run cross validation
        for i_step = 1:n_steps
            i_train = logical(design.train(:,i_step));
            i_test = logical(design.test(:,i_step));
            % pass kernel for current iteration
            kernel_train = kernel(i_train,i_train);
            kernel_test = kernel(i_test,i_train);
            labels_train = design.label(i_train,i_step);
            labels_test = design.label(i_test,i_step);
            % train model with current kernel selection
            model = svmtrain(labels_train,[(1:size(kernel_train,1))',kernel_train],'-s 0 -t 4 -q'); %#ok<*SVMTRAIN>
            % predict all pairwise dvs
            [predicted_labels,~,dv] = svmpredict(labels_test,[(1:size(kernel_test,1))' , kernel_test],model,'-q');

            decoding_out(i_step).predicted_labels = predicted_labels;
            decoding_out(i_step).true_labels = labels_test;
            decoding_out(i_step).decision_values = dv;
            decoding_out(i_step).model = model;
            decoding_out(i_step).opt = [];
        end

        % convert this output to matrix
        out = transres_accuracy_matrix(decoding_out,50);
        result(f,:,:,i_time) = out{1,1};
    end

end

%-average across permutations
result = squeeze(nanmean(result, 1));
result(isnan(result)) = 0;
result = result/100;

rdm = result;

%==========================================================================
else
    
    warning(strcat(mode,' is not implemented yet!'));
    
end

end

%__________________________________________________________________________
function [sigma,shrinkage]=cov1para(x,shrink)

% function sigma=cov1para(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards one-parameter matrix:
%    all variances are the same
%    all covariances are zero
% if shrink is specified, then this value is used for shrinkage

% This version: 04/2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.

% Copyright (c) 2014, Olivier Ledoit and Michael Wolf 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% de-mean returns
[t,n]=size(x);
meanx=mean(x);
x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
sample=(1/t).*(x'*x);

% compute prior
meanvar=mean(diag(sample));
prior=meanvar*eye(n);

if (nargin < 2 | shrink == -1) % compute shrinkage parameters
  
  % what we call p 
  y=x.^2;
  phiMat=y'*y/t-sample.^2;
  phi=sum(sum(phiMat));
  
  % what we call r is not needed for this shrinkage target
  
  % what we call c
  gamma=norm(sample-prior,'fro')^2;

  % compute shrinkage constant
  kappa=phi/gamma;
  shrinkage=max(0,min(1,kappa/t));
    
else % use specified number
  shrinkage=shrink;
end

% compute shrinkage estimator
sigma=shrinkage*prior+(1-shrinkage)*sample;

end

%__________________________________________________________________________
function X = whiten_data(X,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whitening using the Epoch method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_conditions = length(unique(y));
n_sensors = size(X, 2);
n_time = size(X, 3);

sigma_ = nan(n_conditions, n_sensors, n_sensors);
for c = 1:n_conditions
    % compute sigma for each time point, then average across time
    tmp_ = nan(n_time, n_sensors, n_sensors);
    for t = 1:n_time
        tmp_(t, :, :) = cov1para(X(y==c, :, t));
    end
    sigma_(c, :, :) = mean(tmp_, 1);
end
sigma = squeeze(mean(sigma_, 1));  % average across conditions
sigma_inv = sigma^-0.5;

for t = 1:n_time
    X(:, :, t) = squeeze(X(:, :, t)) * sigma_inv;
end

end