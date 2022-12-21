function r_sq_unit = modelVariance(Vc,Vm, D)
% Bharath Talluri, December 2021

% compute rsquared metric between the actual and simulated data

temp_num = (Vc - Vm).^2;
temp_den = (Vc - mean(Vc, 2)).^2;

r_sq_unit.all = 1-(sum(temp_num, 2)./sum(temp_den, 2));
% compute the rsquared for stim and non-stim periods
[~, stim_idx] = intersect(D.timepts2use, D.stim_timepts);
r_sq_unit.stim = 1-(sum(temp_num(:, stim_idx), 2)./sum(temp_den(:, stim_idx), 2));
[~, nonstim_idx] = intersect(D.timepts2use, D.nonstim_timepts);
r_sq_unit.nonstim = 1-(sum(temp_num(:, nonstim_idx), 2)./sum(temp_den(:, nonstim_idx), 2));
% compute the rsquared for inferred and uninferred periods
for i = 1:size(D.robs, 2)
    [~, infer_idx] = intersect(D.timepts2use, find(D.rf_criteria(:, i) == 1));
    r_sq_unit.infer(i, 1) = 1-(sum(temp_num(i, infer_idx), 2)./sum(temp_den(i, infer_idx), 2));
    [~, uninfer_idx] = intersect(D.timepts2use, find(D.rf_criteria(:, i) == 0));
    r_sq_unit.uninfer(i, 1) = 1-(sum(temp_num(i, uninfer_idx), 2)./sum(temp_den(i, uninfer_idx), 2));
    eye_fix_timepts = find(D.rf_criteria(:, i) == 1);
    eye_fix_timepts2use = intersect(D.stim_timepts, eye_fix_timepts);
    [~, fix_idx] = intersect(D.timepts2use, eye_fix_timepts2use);
    r_sq_unit.stim_infer(i, 1) = 1-(sum(temp_num(i, fix_idx), 2)./sum(temp_den(i, fix_idx), 2));
    eye_nonfix_timepts = find(D.rf_criteria(:, i) == 0);
    eye_nonfix_timepts2use = intersect(D.stim_timepts, eye_nonfix_timepts);
    [~, non_fix_idx] = intersect(D.timepts2use, eye_nonfix_timepts2use);
    r_sq_unit.stim_uninfer(i, 1) = 1-(sum(temp_num(i, non_fix_idx), 2)./sum(temp_den(i, non_fix_idx), 2));
    eye_fix_timepts = find(D.rf_criteria(:, i) == 1);
    eye_fix_timepts2use = intersect(D.nonstim_timepts, eye_fix_timepts);
    [~, fix_idx] = intersect(D.timepts2use, eye_fix_timepts2use);
    r_sq_unit.nonstim_infer(i, 1) = 1-(sum(temp_num(i, fix_idx), 2)./sum(temp_den(i, fix_idx), 2));
    eye_nonfix_timepts = find(D.rf_criteria(:, i) == 0);
    eye_nonfix_timepts2use = intersect(D.nonstim_timepts, eye_nonfix_timepts);
    [~, non_fix_idx] = intersect(D.timepts2use, eye_nonfix_timepts2use);
    r_sq_unit.nonstim_uninfer(i, 1) = 1-(sum(temp_num(i, non_fix_idx), 2)./sum(temp_den(i, non_fix_idx), 2));
end
end