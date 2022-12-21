function harvest_vars(varargin)

ip = inputParser();
ip.addParameter('include_eye', 0)
ip.parse(varargin{:});

include_eye = ip.Results.include_eye;
analysis_dir = '../data/analysis';
preprocessed_dir = '../data/preprocessed';
results_dir = '../data/results';
ex_names = {'M1_0008', 'M1_0010', 'M1_0012', 'M1_0013', 'M1_0014', 'M1_0031',...
    'M1_0032', 'M1_0033', 'M1_0034', 'M1_0036', 'M1_0037', 'M1_0038', 'M1_0039',...
    'M1_0040', 'M1_0042', 'M1_0043', 'M1_0044', 'M1_0045', 'M1_0046', 'M1_0048',...
    'M1_0049', 'M1_0050_all1', 'M1_0050_all2', 'M1_0052', 'M1_0054', 'M1_0055',...
    'M1_0056', 'M1_0057', 'M1_0058_15', 'M1_0058_all', 'M1_0059', 'M1_0060',...
    'M1_0063', 'M1_0064', 'M1_0065', 'M1_0066_11', 'M1_0066_all', 'M1_0067_all1',...
    'M1_0067_all2', 'M1_0068', 'M1_0069' 'M1_0071', 'M1_0072_all1', 'M1_0072_all2',...
    'M1_0073_all1', 'M1_0073_all2', 'M1_0074', 'M1_0075', 'M1_0076', 'M1_0077',...
    'M1_0078_all1', 'M1_0078_all2', 'M1_0079', 'M1_0080', 'M2_0803', 'M2_0810',...
    'M2_0822', 'M2_0823', 'M2_0824'};

% define directories
unit_full.all_epochs = [];
unit_task.all_epochs = [];
unit_full.controlled_epochs = [];
unit_task.controlled_epochs = [];
unit_full.uncontrolled_epochs = [];
unit_task.uncontrolled_epochs = [];
unit_full.inferred_epochs = [];
unit_task.inferred_epochs = [];
unit_full.uninferred_epochs = [];
unit_task.uninferred_epochs = [];

face_var_trl_split.rsq.low = [];
face_var_trl_split.rsq.high = [];
body_var_trl_split.rsq.low = [];
body_var_trl_split.rsq.high = [];
face_var_trl_split.mean.low = [];
face_var_trl_split.mean.high = [];
body_var_trl_split.mean.low = [];
body_var_trl_split.mean.high = [];
face_var_trl_split.var.low = [];
face_var_trl_split.var.high = [];
body_var_trl_split.var.low = [];
body_var_trl_split.var.high = [];

all_epochs.drift_unique = [];
all_epochs.move_unique = [];
all_epochs.body_unique = [];
all_epochs.face_unique = [];
all_epochs.task_unique = [];
all_epochs.pupil_unique = [];
all_epochs.eye_unique = [];

controlled_epochs.drift_unique = [];
controlled_epochs.move_unique = [];
controlled_epochs.body_unique = [];
controlled_epochs.face_unique = [];
controlled_epochs.task_unique = [];
controlled_epochs.pupil_unique = [];
controlled_epochs.eye_unique = [];

uncontrolled_epochs.drift_unique = [];
uncontrolled_epochs.move_unique = [];
uncontrolled_epochs.body_unique = [];
uncontrolled_epochs.face_unique = [];
uncontrolled_epochs.task_unique = [];
uncontrolled_epochs.pupil_unique = [];
uncontrolled_epochs.eye_unique = [];

uncontrolled_epochs.infer.drift_unique = [];
uncontrolled_epochs.infer.move_unique = [];
uncontrolled_epochs.infer.body_unique = [];
uncontrolled_epochs.infer.face_unique = [];
uncontrolled_epochs.infer.task_unique = [];
uncontrolled_epochs.infer.pupil_unique = [];
uncontrolled_epochs.infer.eye_unique = [];

uncontrolled_epochs.uninfer.drift_unique = [];
uncontrolled_epochs.uninfer.move_unique = [];
uncontrolled_epochs.uninfer.body_unique = [];
uncontrolled_epochs.uninfer.face_unique = [];
uncontrolled_epochs.uninfer.task_unique = [];
uncontrolled_epochs.uninfer.pupil_unique = [];
uncontrolled_epochs.uninfer.eye_unique = [];

inferred_epochs.drift_unique = [];
inferred_epochs.move_unique = [];
inferred_epochs.body_unique = [];
inferred_epochs.face_unique = [];
inferred_epochs.task_unique = [];
inferred_epochs.pupil_unique = [];
inferred_epochs.eye_unique = [];

uninferred_epochs.drift_unique = [];
uninferred_epochs.move_unique = [];
uninferred_epochs.body_unique = [];
uninferred_epochs.face_unique = [];
uninferred_epochs.task_unique = [];
uninferred_epochs.pupil_unique = [];
uninferred_epochs.eye_unique = [];

group_mean_rates = [];
group_mean_quartile_rates = [];
group_rf_criteria = [];
group_rf_inferred = [];
group_rf_uninferred = [];
group_sv.body_var = [];
group_sv.face_var = [];
group_sv.body_mean = [];
group_sv.face_mean = [];
exp_idx = [];
sess_units = [];
psth_rsq = [];
for ex = 1:length(ex_names)
    exname = ex_names{ex};
    if ~include_eye
        load(sprintf('%s/%s/preprocessed_data.mat', preprocessed_dir, exname), 'D');
    else
        load(sprintf('%s/%s/preprocessed_data.mat', preprocessed_dir, exname), 'D_eye');
        D = D_eye;clear D_eye;
    end
    num_units = zeros(length(D.units), 1);
    if isfield(D, 'area')
        if isfield(D.area, 'V1_idx') && ~isempty(D.area.V1_idx)
            num_units(D.area.V1_idx) = 1;
        end
        if isfield(D.area, 'V2_idx') && ~isempty(D.area.V2_idx)
            num_units(D.area.V2_idx) = 2;
        end
        if isfield(D.area, 'V3_idx') && ~isempty(D.area.V3_idx)
            num_units(D.area.V3_idx) = 3;
        end
    end
    sess_units = [sess_units;num_units];
    % load the analysis variables
    if ~ include_eye
        load(sprintf('%s/%s/analysisVars.mat', analysis_dir, exname), 'Vc', 'regIdx', 'regLabels');
    else
        load(sprintf('%s/%s/analysisVars_withEye.mat', analysis_dir, exname), 'Vc', 'regIdx', 'regLabels');
    end

    % get face and body SV information
    num_sv = length(find(regIdx == find(strcmp(regLabels, 'face_svd'))));
    vidR = [D.faceSVD(:, 1:num_sv) D.bodySVD(:, 1:num_sv)];
    % remove repeated/invalid timepoints from the analysis
    vidR(D.timepts2remove, :) = [];
    % zscore video SV
    vidR = zscore(double(vidR));
    trial_start = D.blocks_clean(:, 1);
    trial_end = D.blocks_clean(:, 2) - 1;
    face_sv_var = nan(num_sv, 200);
    face_sv_mean = nan(num_sv, 200);
    for i = 1:num_sv
        face_sv_trl = nan(length(trial_start), 200);
        for j = 1:length(trial_start)
            face_sv_trl(j, 1:min(trial_end(j)-trial_start(j)+1, 200)) = vidR(trial_start(j):min(trial_end(j), trial_start(j)+199), i)';
        end
        face_sv_var(i, :) = nanvar(face_sv_trl);
        face_sv_mean(i, :) = nanmean(abs(face_sv_trl));
    end
    body_sv_var = nan(num_sv, 200);
    body_sv_mean = nan(num_sv, 200);
    for i = 1:num_sv
        body_sv_trl = nan(length(trial_start), 200);
        for j = 1:length(trial_start)
            body_sv_trl(j, 1:min(trial_end(j)-trial_start(j)+1, 200)) = vidR(trial_start(j):min(trial_end(j), trial_start(j)+199), i+num_sv)';
        end
        body_sv_var(i, :) = nanvar(body_sv_trl);
        body_sv_mean(i, :) = nanmean(abs(body_sv_trl));
    end

    group_sv.body_var = [group_sv.body_var; nanmean(body_sv_var, 1)];
    group_sv.face_var = [group_sv.face_var; nanmean(face_sv_var, 1)];
    group_sv.body_mean = [group_sv.body_mean; nanmean(body_sv_mean, 1)];
    group_sv.face_mean = [group_sv.face_mean; nanmean(face_sv_mean, 1)];

    % get the mean stimulus-driven activity per unit- we will use this
    % information later to select 'dead' units
    [mean_rate, mean_quartile_rate] = unit_activity(Vc, D);
    group_mean_rates = [group_mean_rates; mean_rate];
    group_mean_quartile_rates = [group_mean_quartile_rates; mean_quartile_rate];
    group_rf_criteria = [group_rf_criteria; sum(D.rf_criteria, 1)'];
    % get the proportion of inferred and uninferred retinal input
    % timepoints
    for i = 1:size(D.robs, 2)
        eye_fix_timepts = find(D.rf_criteria(:, i) == 1);
        [~, fix_idx] = intersect(D.timepts2use, eye_fix_timepts);
        group_rf_inferred = [group_rf_inferred; length(fix_idx)];

        eye_nonfix_timepts = find(D.rf_criteria(:, i) == 0);
        [~, non_fix_idx] = intersect(D.timepts2use, eye_nonfix_timepts);
        group_rf_uninferred = [group_rf_uninferred; length(non_fix_idx)];
    end

    if ~include_eye
        load(sprintf('%s/%s/crossval_modelFits.mat', analysis_dir, exname), 'Vfull', 'VtaskMod', 'fullVar', 'taskModVar');
    else
        load(sprintf('%s/%s/crossval_modelFits_withEye.mat', analysis_dir, exname), 'Vfull', 'VtaskMod', 'fullVar', 'taskModVar');
    end
    % get the R2 for unit mean and prediucted a ctivity by full model
    spk_psth = spk_timeseries(Vc, D.blocks_clean);
    fullModel_psth = spk_timeseries(Vfull, D.blocks_clean);
    for this_idx = 1:length(D.units)
        r = this_rsq(spk_psth(this_idx, :)', fullModel_psth(this_idx, :)');
        psth_rsq = [psth_rsq; r];
    end

    % compute the variance across SVs per trial during controlled epochs
    [~, stim_idx] = intersect(D.timepts2use, D.stim_timepts);
    face_var_controlled = nan(length(trial_start), num_sv);
    body_var_controlled = nan(length(trial_start), num_sv);
    face_mean_controlled = nan(length(trial_start), num_sv);
    body_mean_controlled = nan(length(trial_start), num_sv);
    for i = 1:num_sv
        for j = 1:length(trial_start)
            controlled_trial_idx = intersect(trial_start(j):trial_end(j), stim_idx);
            face_var_controlled(j, i) = nanvar(vidR(controlled_trial_idx, i));
            body_var_controlled(j, i) = nanvar(vidR(controlled_trial_idx, i+num_sv));
            face_mean_controlled(j, i) = nanmean(abs(vidR(controlled_trial_idx, i)));
            body_mean_controlled(j, i) = nanmean(abs(vidR(controlled_trial_idx, i+num_sv)));
        end
    end
    trl_face_var_controlled = nanmean(face_var_controlled, 2);
    trl_body_var_controlled = nanmean(body_var_controlled, 2);
    trl_face_mean_controlled = nanmean(face_mean_controlled, 2);
    trl_body_mean_controlled = nanmean(body_mean_controlled, 2);

    high_face_var_trls = find(trl_face_var_controlled >= prctile(trl_face_var_controlled, 200/3));
    low_face_var_trls = find(trl_face_var_controlled <= prctile(trl_face_var_controlled, 100/3));
    high_face_idx = [];
    for j = high_face_var_trls'
        controlled_trial_idx = intersect(trial_start(j):trial_end(j), stim_idx);
        high_face_idx = [high_face_idx; controlled_trial_idx];
    end
    low_face_idx = [];
    for j = low_face_var_trls'
        controlled_trial_idx = intersect(trial_start(j):trial_end(j), stim_idx);
        low_face_idx = [low_face_idx; controlled_trial_idx];
    end
    for this_idx = 1:length(D.units)
        face_var_trl_split.rsq.high = [face_var_trl_split.rsq.high; this_rsq(Vc(this_idx, high_face_idx), Vfull(this_idx, high_face_idx))-this_rsq(Vc(this_idx, high_face_idx), VtaskMod(this_idx, high_face_idx))];
        face_var_trl_split.rsq.low = [face_var_trl_split.rsq.low;this_rsq(Vc(this_idx, low_face_idx), Vfull(this_idx, low_face_idx))-this_rsq(Vc(this_idx, low_face_idx), VtaskMod(this_idx, low_face_idx))];
    end
    face_var_trl_split.mean.high = [face_var_trl_split.mean.high; mean(trl_face_mean_controlled(high_face_var_trls))];
    face_var_trl_split.mean.low = [face_var_trl_split.mean.low;mean(trl_face_mean_controlled(low_face_var_trls))];
    face_var_trl_split.var.high = [face_var_trl_split.var.high; mean(trl_face_var_controlled(high_face_var_trls))];
    face_var_trl_split.var.low = [face_var_trl_split.var.low;mean(trl_face_var_controlled(low_face_var_trls))];


    high_body_var_trls = find(trl_body_var_controlled >= prctile(trl_body_var_controlled, 200/3));
    low_body_var_trls = find(trl_body_var_controlled <= prctile(trl_body_var_controlled, 100/3));
    high_body_idx = [];
    for j = high_body_var_trls'
        controlled_trial_idx = intersect(trial_start(j):trial_end(j), stim_idx);
        high_body_idx = [high_body_idx; controlled_trial_idx];
    end
    low_body_idx = [];
    for j = low_body_var_trls'
        controlled_trial_idx = intersect(trial_start(j):trial_end(j), stim_idx);
        low_body_idx = [low_body_idx; controlled_trial_idx];
    end
    for this_idx = 1:length(D.units)
        body_var_trl_split.rsq.high = [body_var_trl_split.rsq.high;this_rsq(Vc(this_idx, high_body_idx), Vfull(this_idx, high_body_idx))-this_rsq(Vc(this_idx, high_body_idx), VtaskMod(this_idx, high_body_idx))];
        body_var_trl_split.rsq.low = [body_var_trl_split.rsq.low;this_rsq(Vc(this_idx, low_body_idx), Vfull(this_idx, low_body_idx))-this_rsq(Vc(this_idx, low_body_idx), VtaskMod(this_idx, low_body_idx))];
    end
    body_var_trl_split.mean.high = [body_var_trl_split.mean.high; mean(trl_body_mean_controlled(high_body_var_trls))];
    body_var_trl_split.mean.low = [body_var_trl_split.mean.low;mean(trl_body_mean_controlled(low_body_var_trls))];
    body_var_trl_split.var.high = [body_var_trl_split.var.high; mean(trl_body_var_controlled(high_body_var_trls))];
    body_var_trl_split.var.low = [body_var_trl_split.var.low;mean(trl_body_var_controlled(low_body_var_trls))];

    if include_eye
        load(sprintf('%s/%s/uniqueVariance_withEye.mat', analysis_dir, exname), 'drift_unique', 'body_unique', 'face_unique', 'eye_unique', 'pupil_unique', 'task_unique');
    else
        load(sprintf('%s/%s/uniqueVariance.mat', analysis_dir, exname), 'drift_unique', 'body_unique', 'face_unique', 'task_unique');
    end

    unit_full.all_epochs = [unit_full.all_epochs; fullVar.all];
    unit_task.all_epochs = [unit_task.all_epochs; taskModVar.all];
    unit_full.controlled_epochs = [unit_full.controlled_epochs; fullVar.stim];
    unit_task.controlled_epochs = [unit_task.controlled_epochs; taskModVar.stim];
    unit_full.uncontrolled_epochs = [unit_full.uncontrolled_epochs; fullVar.nonstim];
    unit_task.uncontrolled_epochs = [unit_task.uncontrolled_epochs; taskModVar.nonstim];
    unit_full.inferred_epochs = [unit_full.inferred_epochs; fullVar.infer];
    unit_task.inferred_epochs = [unit_task.inferred_epochs; taskModVar.infer];
    unit_full.uninferred_epochs = [unit_full.uninferred_epochs; fullVar.uninfer];
    unit_task.uninferred_epochs = [unit_task.uninferred_epochs; taskModVar.uninfer];

    all_epochs.drift_unique = [all_epochs.drift_unique; drift_unique.all];
    all_epochs.move_unique = [all_epochs.move_unique; fullVar.all - taskModVar.all];
    all_epochs.body_unique = [all_epochs.body_unique; body_unique.all];
    all_epochs.face_unique = [all_epochs.face_unique; face_unique.all];
    all_epochs.task_unique = [all_epochs.task_unique; task_unique.all];
    if include_eye
        all_epochs.pupil_unique = [all_epochs.pupil_unique; pupil_unique.all];
        all_epochs.eye_unique = [all_epochs.eye_unique; eye_unique.all];
    end

    controlled_epochs.drift_unique = [controlled_epochs.drift_unique; drift_unique.stim];
    controlled_epochs.move_unique = [controlled_epochs.move_unique; fullVar.stim - taskModVar.stim];
    controlled_epochs.body_unique = [controlled_epochs.body_unique; body_unique.stim];
    controlled_epochs.face_unique = [controlled_epochs.face_unique; face_unique.stim];
    controlled_epochs.task_unique = [controlled_epochs.task_unique; task_unique.stim];
    if include_eye
        controlled_epochs.pupil_unique = [controlled_epochs.pupil_unique; pupil_unique.stim];
        controlled_epochs.eye_unique = [controlled_epochs.eye_unique; eye_unique.stim];
    end

    uncontrolled_epochs.drift_unique = [uncontrolled_epochs.drift_unique; drift_unique.nonstim];
    uncontrolled_epochs.move_unique = [uncontrolled_epochs.move_unique; fullVar.nonstim - taskModVar.nonstim];
    uncontrolled_epochs.body_unique = [uncontrolled_epochs.body_unique; body_unique.nonstim];
    uncontrolled_epochs.face_unique = [uncontrolled_epochs.face_unique; face_unique.nonstim];
    uncontrolled_epochs.task_unique = [uncontrolled_epochs.task_unique; task_unique.nonstim];
    if include_eye
        uncontrolled_epochs.pupil_unique = [uncontrolled_epochs.pupil_unique; pupil_unique.nonstim];
        uncontrolled_epochs.eye_unique = [uncontrolled_epochs.eye_unique; eye_unique.nonstim];
    end

    inferred_epochs.drift_unique = [inferred_epochs.drift_unique; drift_unique.infer];
    inferred_epochs.move_unique = [inferred_epochs.move_unique; fullVar.infer - taskModVar.infer];
    inferred_epochs.body_unique = [inferred_epochs.body_unique; body_unique.infer];
    inferred_epochs.face_unique = [inferred_epochs.face_unique; face_unique.infer];
    inferred_epochs.task_unique = [inferred_epochs.task_unique; task_unique.infer];
    if include_eye
        inferred_epochs.pupil_unique = [inferred_epochs.pupil_unique; pupil_unique.infer];
        inferred_epochs.eye_unique = [inferred_epochs.eye_unique; eye_unique.infer];
    end

    uninferred_epochs.drift_unique = [uninferred_epochs.drift_unique; drift_unique.uninfer];
    uninferred_epochs.move_unique = [uninferred_epochs.move_unique; fullVar.uninfer - taskModVar.uninfer];
    uninferred_epochs.body_unique = [uninferred_epochs.body_unique; body_unique.uninfer];
    uninferred_epochs.face_unique = [uninferred_epochs.face_unique; face_unique.uninfer];
    uninferred_epochs.task_unique = [uninferred_epochs.task_unique; task_unique.uninfer];
    if include_eye
        uninferred_epochs.pupil_unique = [uninferred_epochs.pupil_unique; pupil_unique.uninfer];
        uninferred_epochs.eye_unique = [uninferred_epochs.eye_unique; eye_unique.uninfer];
    end

    uncontrolled_epochs.infer.drift_unique = [uncontrolled_epochs.infer.drift_unique; drift_unique.nonstim_infer];
    uncontrolled_epochs.infer.move_unique = [uncontrolled_epochs.infer.move_unique; fullVar.nonstim_infer - taskModVar.nonstim_infer];
    uncontrolled_epochs.infer.body_unique = [uncontrolled_epochs.infer.body_unique; body_unique.nonstim_infer];
    uncontrolled_epochs.infer.face_unique = [uncontrolled_epochs.infer.face_unique; face_unique.nonstim_infer];
    uncontrolled_epochs.infer.task_unique = [uncontrolled_epochs.infer.task_unique; task_unique.nonstim_infer];
    if include_eye
        uncontrolled_epochs.infer.pupil_unique = [uncontrolled_epochs.infer.pupil_unique; pupil_unique.nonstim_infer];
        uncontrolled_epochs.infer.eye_unique = [uncontrolled_epochs.infer.eye_unique; eye_unique.nonstim_infer];
    end

    uncontrolled_epochs.uninfer.drift_unique = [uncontrolled_epochs.uninfer.drift_unique; drift_unique.nonstim_uninfer];
    uncontrolled_epochs.uninfer.move_unique = [uncontrolled_epochs.uninfer.move_unique; fullVar.nonstim_uninfer - taskModVar.nonstim_uninfer];
    uncontrolled_epochs.uninfer.body_unique = [uncontrolled_epochs.uninfer.body_unique; body_unique.nonstim_uninfer];
    uncontrolled_epochs.uninfer.face_unique = [uncontrolled_epochs.uninfer.face_unique; face_unique.nonstim_uninfer];
    uncontrolled_epochs.uninfer.task_unique = [uncontrolled_epochs.uninfer.task_unique; task_unique.nonstim_uninfer];
    if include_eye
        uncontrolled_epochs.uninfer.pupil_unique = [uncontrolled_epochs.uninfer.pupil_unique; pupil_unique.nonstim_uninfer];
        uncontrolled_epochs.uninfer.eye_unique = [uncontrolled_epochs.uninfer.eye_unique; eye_unique.nonstim_uninfer];
    end

    exp_idx = [exp_idx; ex*ones(length(num_units), 1)];
end
if include_eye
    save(sprintf('%s/harvested_vars_withEye.mat', results_dir), 'group_sv', 'group_mean_rates', 'group_rf_criteria', 'group_rf_uninferred', 'group_rf_inferred', 'group_mean_quartile_rates', 'body_var_trl_split', 'face_var_trl_split', 'psth_rsq', 'unit_full', 'unit_task', 'sess_units', 'all_epochs', 'controlled_epochs', 'uncontrolled_epochs', 'inferred_epochs', 'uninferred_epochs', 'exp_idx', 'ex_names', '-v7.3');
else
    save(sprintf('%s/harvested_vars.mat', results_dir), 'group_sv', 'group_mean_rates', 'group_rf_criteria', 'group_rf_uninferred', 'group_rf_inferred', 'group_mean_quartile_rates', 'body_var_trl_split', 'face_var_trl_split', 'psth_rsq', 'unit_full', 'unit_task', 'sess_units', 'all_epochs', 'controlled_epochs', 'uncontrolled_epochs', 'inferred_epochs', 'uninferred_epochs', 'exp_idx', 'ex_names', '-v7.3');
end
end

function Vt = trial_spikes(Vs, blocks)
% get the trial-related spiking activity
Vt = nan(size(blocks, 1), max(blocks(:, 2) - blocks(:, 1)));
num_trials = size(blocks, 1);
for tr = 1:num_trials
    Vt(tr, 1:blocks(tr,2)-blocks(tr,1)) = Vs(blocks(tr,1):blocks(tr,2)-1);
end
end

function spks = spk_timeseries(Vc, blocks)
% get spike densitiy functions
spks = cell(size(Vc, 1), 1);
kwidth = 4.8; % SD of smoothing kernel width in ms
x = 0:kwidth*3;
half_g = exp(-(x.^2)/(2 * kwidth^2));
scale = 60/sum(half_g); % ms resolution for spikes
for c = 1:size(Vc, 1) % unit
    unit_spk = Vc(c, :);
    unit_spk_trl = trial_spikes(unit_spk, blocks);
    spks{c, 1} = conv(mean(unit_spk_trl, 1),half_g,'same')*scale;
end
spks = cell2mat(spks(:, 1));
end

function r_sq_unit = this_rsq(Vc,Vm)
% compute rsquared metric between the actual and simulated data
fin_idx = find(~isnan(Vc) & ~isnan(Vm));
temp_num = (Vc(fin_idx) - Vm(fin_idx)).^2;
temp_den = (Vc(fin_idx) - mean(Vc(fin_idx))).^2;
r_sq_unit = 1-(sum(temp_num)./sum(temp_den));
end