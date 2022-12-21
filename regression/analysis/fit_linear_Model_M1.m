function fit_linear_Model_M1(exname, include_eye)

results_dir = sprintf('../data/analysis/%s', exname);

% load preprocessed dataset
if ~include_eye
    load(sprintf('../data/preprocessed/%s/preprocessed_data.mat', exname), 'D');
else
    load(sprintf('../data/preprocessed/%s/preprocessed_data.mat', exname), 'D_eye');
    D = D_eye;clear D_eye;
end

if ~include_eye
    load(sprintf('%s/analysisVars.mat', results_dir), 'Vc', 'opts', 'fullR', 'regIdx', 'regLabels', 'ModR', 'dataIdx');
else
    load(sprintf('%s/analysisVars_withEye.mat', results_dir), 'Vc', 'opts', 'fullR', 'regIdx', 'regLabels', 'ModR', 'dataIdx');
end
%% fit model to spike data
disp("Running Ridge Regression")
[Vridge, ridgeBeta, ~, ~, ~, ~] = crossValModel(fullR, Vc, regLabels, regIdx, regLabels, 1);
ridgeBeta = cell2mat(ridgeBeta);
disp("Done")

%% Compute variance explained for the full model
fullVar_Ridge = modelVariance(Vc,Vridge,D); %compute explained variance

%% run cross-validation

disp("Running Cross validation")

%full model
[Vfull, full_crossValBeta, full_crossvalR, full_crossvalIdx, ~, fullLabels] = crossValModel(fullR, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);
fullVar = modelVariance(Vc,Vfull,D);  %compute explained variance

% get the drift variable prediction
targetGroups = {'drift'};
cvars = find(contains(fullLabels, targetGroups));
cIdx = full_crossvalIdx==cvars;
cIdx = any(cIdx, 2);

Vfull_drift = zeros(size(Vc),'single'); % pre-allocate
for iFolds = 1:opts.folds
    Vfull_drift(:,~dataIdx(iFolds, :)) = ((full_crossvalR(~dataIdx(iFolds, :),cIdx) - mean(full_crossvalR(~dataIdx(iFolds, :),cIdx), 1)) * full_crossValBeta{iFolds}(cIdx, :))';
end

% task model with modulation
[VtaskMod, taskMod_crossValBeta, ~, ~, ~, ~] = crossValModel(ModR, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);
taskModVar = modelVariance(Vc,VtaskMod,D);

disp("Done")

if ~include_eye
save(sprintf('%s/crossval_modelFits.mat', results_dir), 'fullVar_Ridge', 'ridgeBeta', 'Vridge', 'dataIdx', 'Vfull', 'Vfull_drift', 'VtaskMod', 'full_crossValBeta', 'taskMod_crossValBeta', 'fullVar', 'taskModVar', '-v7.3');
else
save(sprintf('%s/crossval_modelFits_withEye.mat', results_dir), 'fullVar_Ridge', 'ridgeBeta', 'Vridge', 'dataIdx', 'Vfull', 'Vfull_drift', 'VtaskMod', 'full_crossValBeta', 'taskMod_crossValBeta', 'fullVar', 'taskModVar', '-v7.3');
end
%% compute unique variance metrics

disp("Computing unique variance")

% drift regressors
shuffle_regressors = 'drift';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
[V_driftUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);

% body regressors
shuffle_regressors = 'body_svd';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
[V_bodyUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);

% face regressors
shuffle_regressors = 'face_svd';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
[V_faceUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);

% task regressors
svd_idx = find(contains(fullLabels, {'svd', 'eye', 'pupil', 'drift'}));
taskLabels = setdiff(fullLabels, fullLabels(svd_idx));
task_idx = find(contains(fullLabels, taskLabels));
use_labels = fullLabels(task_idx);
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, use_labels))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
[V_taskUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);

if include_eye
    % eye regressors
    shuffle_regressors = 'eye';
    shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
    shuffled_R = fullR;
    for idx = shuffle_idx'
        shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
    end
    [V_eyeUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);

    % pupil regressors
    shuffle_regressors = 'pupil';
    shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
    shuffled_R = fullR;
    for idx = shuffle_idx'
        shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
    end
    [V_pupilUnique, ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc, regLabels, regIdx, regLabels, opts.folds, dataIdx);
end
disp("Done")

% compute the unique variance explained metrics
fprintf("Compute variance explained for %s \n", exname);
shuffle_task = modelVariance(Vc, V_taskUnique, D);
f = fieldnames(fullVar);
for k=1:numel(f)
    task_unique.(f{k}) = fullVar.(f{k})-shuffle_task.(f{k});
end
shuffle_drift = modelVariance(Vc, V_driftUnique, D);
f = fieldnames(fullVar);
for k=1:numel(f)
    drift_unique.(f{k}) = fullVar.(f{k})-shuffle_drift.(f{k});
end
shuffle_body = modelVariance(Vc, V_bodyUnique, D);
f = fieldnames(fullVar);
for k=1:numel(f)
    body_unique.(f{k}) = fullVar.(f{k})-shuffle_body.(f{k});
end
shuffle_face = modelVariance(Vc, V_faceUnique, D);
f = fieldnames(fullVar);
for k=1:numel(f)
    face_unique.(f{k}) = fullVar.(f{k})-shuffle_face.(f{k});
end

if include_eye
    shuffle_pupil = modelVariance(Vc, V_pupilUnique, D);
    f = fieldnames(fullVar);
    for k=1:numel(f)
        pupil_unique.(f{k}) = fullVar.(f{k})-shuffle_pupil.(f{k});
    end
    shuffle_eye = modelVariance(Vc, V_eyeUnique, D);
    f = fieldnames(fullVar);
    for k=1:numel(f)
        eye_unique.(f{k}) = fullVar.(f{k})-shuffle_eye.(f{k});
    end
end
disp("Done")

% save results
if include_eye
    save(sprintf('%s/uniqueVariance_withEye.mat', results_dir), 'drift_unique', 'body_unique', 'face_unique', 'eye_unique', 'pupil_unique', 'task_unique', '-v7.3');
else
    save(sprintf('%s/uniqueVariance.mat', results_dir), 'drift_unique', 'body_unique', 'face_unique', 'task_unique', '-v7.3');
end