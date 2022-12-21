function fit_linear_Model_M2(exname, include_eye)

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
% for animal M2, we will fit the models taking the attention, and
% receptive field location into account

disp("Running Cross validation")

Vfull = nan(size(Vc));
VtaskMod = nan(size(Vc));
Vridge = nan(size(Vc));
Vfull_drift = nan(size(Vc));
ridgeBeta = nan(size(fullR, 2), size(Vc, 1));

% fit separately for units with rfs on contralateral hemifield, and remove
% regressors modelling stimulus on the ipsilateral side
rf_locs = {'left', 'right'};

for rf_loc = rf_locs
    rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
    if isempty (rf_idx)
        continue
    else
        nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
        RF_labels = setdiff(regLabels, nonRF_labels);
        [Vfull(rf_idx, :), full_crossValBeta, full_crossvalR, full_crossvalIdx, ~, fullLabels] = crossValModel(fullR, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
        [VtaskMod(rf_idx, :), taskMod_crossValBeta, ~, ~, ~, ~] = crossValModel(ModR, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
        [Vridge(rf_idx, :), ridge_beta, ~, ~, ~, ~] = crossValModel(fullR, Vc(rf_idx, :), regLabels, regIdx, regLabels, 1);
        ridgeBeta(:, rf_idx) = cell2mat(ridge_beta);
        % get the drift variable prediction
        targetGroups = {'drift'};
        cvars = find(contains(fullLabels, targetGroups));
        cIdx = full_crossvalIdx==cvars;
        cIdx = any(cIdx, 2);
        for iFolds = 1:opts.folds
            Vfull_drift(rf_idx,~dataIdx(iFolds, :)) = ((full_crossvalR(~dataIdx(iFolds, :),cIdx) - mean(full_crossvalR(~dataIdx(iFolds, :),cIdx), 1)) * full_crossValBeta{iFolds}(cIdx, :))';
        end
    end
end
%% Compute variance explained for the full model
fullVar_Ridge = modelVariance(Vc,Vridge,D); %compute explained variance
fullVar = modelVariance(Vc,Vfull,D);
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
for rf_loc = rf_locs
    rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
    if isempty (rf_idx)
        continue
    else
        nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
        RF_labels = setdiff(regLabels, nonRF_labels);
        [V_driftUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
    end
end

% body regressors
shuffle_regressors = 'body_svd';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
for rf_loc = rf_locs
    rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
    if isempty (rf_idx)
        continue
    else
        nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
        RF_labels = setdiff(regLabels, nonRF_labels);
        [V_bodyUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
    end
end

% face regressors
shuffle_regressors = 'face_svd';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
shuffled_R = fullR;
for idx = shuffle_idx'
    shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
end
for rf_loc = rf_locs
    rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
    if isempty (rf_idx)
        continue
    else
        nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
        RF_labels = setdiff(regLabels, nonRF_labels);
        [V_faceUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
    end
end

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
for rf_loc = rf_locs
    rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
    if isempty (rf_idx)
        continue
    else
        nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
        RF_labels = setdiff(regLabels, nonRF_labels);
        [V_taskUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
    end
end

if include_eye
    % eye regressors
    shuffle_regressors = 'eye';
    shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
    shuffled_R = fullR;
    for idx = shuffle_idx'
        shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
    end
    for rf_loc = rf_locs
        rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
        if isempty (rf_idx)
            continue
        else
            nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
            RF_labels = setdiff(regLabels, nonRF_labels);
            [V_eyeUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
        end
    end

    % pupil regressors
    shuffle_regressors = 'pupil';
    shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
    shuffled_R = fullR;
    for idx = shuffle_idx'
        shuffled_R(:, idx) = datasample(shuffled_R(:, idx), size(shuffled_R, 1), 1);
    end
    for rf_loc = rf_locs
        rf_idx = find(strcmp(rf_loc, {D.units.rf_loc}));
        if isempty (rf_idx)
            continue
        else
            nonRF_labels = regLabels(find(contains(regLabels, setdiff(rf_locs, rf_loc))));
            RF_labels = setdiff(regLabels, nonRF_labels);
            [V_pupilUnique(rf_idx, :), ~, ~, ~, ~, ~] = crossValModel(shuffled_R, Vc(rf_idx, :), RF_labels, regIdx, regLabels, opts.folds, dataIdx);
        end
    end
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