function fig2(animal, varargin)

ip = inputParser();
ip.addParameter('num_sv', 30)
ip.addParameter('time_pre_stim', .3)
ip.addParameter('time_per_trial', 3)
ip.addParameter('dt', 1/60)
ip.addParameter('include_eye', 0)
ip.parse(varargin{:});

include_eye = ip.Results.include_eye;
num_sv = ip.Results.num_sv;
dt = ip.Results.dt;
time_per_trial = ip.Results.time_per_trial;
time_pre_stim = ip.Results.time_pre_stim;

% define directories
results_dir = '../data/results';
if ~include_eye
    load(sprintf('%s/harvested_vars.mat', results_dir));
else
    load(sprintf('%s/harvested_vars_withEye.mat', results_dir));
end

figure(1);clf;
colormap linspecer;
cols = linspecer(9,'qualitative');
cols(1, :) = [0 0 1];
cols(2, :) = [0 0.65 1];

if strcmp(animal, 'M2')
    good_units = find(unit_full.controlled_epochs >= 0 & unit_full.uncontrolled_epochs >= 0 & all(group_mean_quartile_rates >= 2, 2) & group_mean_rates >= 2 & exp_idx > 54); %
elseif strcmp(animal, 'M1')
    good_units = find(unit_full.controlled_epochs >= 0 & unit_full.uncontrolled_epochs >= 0 & all(group_mean_quartile_rates >= 2, 2) & group_mean_rates >= 2 & exp_idx < 55); %
else
    good_units = find(unit_full.controlled_epochs >= 0 & unit_full.uncontrolled_epochs >= 0 & all(group_mean_quartile_rates >= 2, 2) & group_mean_rates >= 2); %
    load(sprintf('%s/selected_units.mat', results_dir));
end

%% plot SV mean and variance as a function of trial time, fig 2B

trial_time = -time_pre_stim+dt:dt:time_per_trial;
x = trial_time;
idx = 1:2;
x(idx) = [];

% animal M1
m = nanmean(group_sv.body_var(1:54, 1:198),1);
sd = nanstd(group_sv.body_var(1:54, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_body = x; m_body = m;sd_body = sd;
m = nanmean(group_sv.face_var(1:54, 1:198),1);
sd = nanstd(group_sv.face_var(1:54, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_face = x; m_face = m;sd_face = sd;
subplot(8,4,1);hold on;
myshadedErrorBar(x_body, m_body, sd_body, '-', cols(1, :), 0.4);
myshadedErrorBar(x_face, m_face, sd_face, '-', cols(2, :), 0.4);
h1 = line(x_body,m_body,'Color',cols(1, :),'LineWidth',2);
h2 = line(x_face,m_face,'Color',cols(2, :),'LineWidth',2);
vertline(0, 'k--');legend([h1, h2], 'Body SVD', 'Face SVD', 'Location', 'Best', 'FontSize', 6);legend('boxoff');
ylabel('Variance', 'FontSize', 6);title('M1', 'FontSize', 6)
set(gca, 'XLim', [-time_pre_stim, time_per_trial], 'XTick', [-time_pre_stim 0:0.5:time_per_trial], 'YLim', [0 4], 'YTick', 0:4, 'YTickLabel', [], 'FontSize', 5);offsetAxes;

m = nanmean(group_sv.body_mean(1:54, 1:198),1);
sd = nanstd(group_sv.body_mean(1:54, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_body = x; m_body = m;sd_body = sd;
m = nanmean(group_sv.face_mean(1:54, 1:198),1);
sd = nanstd(group_sv.face_mean(1:54, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_face = x; m_face = m;sd_face = sd;
subplot(8,4,2);hold on;
myshadedErrorBar(x_body, m_body, sd_body, '-', cols(1, :), 0.4);
myshadedErrorBar(x_face, m_face, sd_face, '-', cols(2, :), 0.4);
line(x_body,m_body,'Color',cols(1, :),'LineWidth',2);
line(x_face,m_face,'Color',cols(2, :),'LineWidth',2);
vertline(0, 'k--');
ylabel('Mean', 'FontSize', 6);title('M1', 'FontSize', 6)
set(gca, 'XLim', [-time_pre_stim, time_per_trial], 'XTick', [-time_pre_stim 0:0.5:time_per_trial], 'YLim', [0 4], 'YTick', 0:4, 'YTickLabel', [], 'FontSize', 5);offsetAxes;

% animal M2
m = nanmean(group_sv.body_var(55:end, 1:198),1);
sd = nanstd(group_sv.body_var(55:end, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_body = x; m_body = m;sd_body = sd;
m = nanmean(group_sv.face_var(55:end, 1:198),1);
sd = nanstd(group_sv.face_var(55:end, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_face = x; m_face = m;sd_face = sd;
subplot(8,4,5);hold on;
myshadedErrorBar(x_body, m_body, sd_body, '-', cols(1, :), 0.4);
myshadedErrorBar(x_face, m_face, sd_face, '-', cols(2, :), 0.4);
line(x_body,m_body,'Color',cols(1, :),'LineWidth',2);
line(x_face,m_face,'Color',cols(2, :),'LineWidth',2);
vertline(0, 'k--');
ylabel('Variance', 'FontSize', 6);title('M2', 'FontSize', 6)
set(gca, 'XLim', [-time_pre_stim, time_per_trial], 'XTick', [-time_pre_stim 0:0.5:time_per_trial], 'YLim', [0 4], 'YTick', 0:4, 'YTickLabel', [], 'FontSize', 5);offsetAxes;

m = nanmean(group_sv.body_mean(55:end, 1:198),1);
sd = nanstd(group_sv.body_mean(55:end, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_body = x; m_body = m;sd_body = sd;
m = nanmean(group_sv.face_mean(55:end, 1:198),1);
sd = nanstd(group_sv.face_mean(55:end, 1:198),0,1);
m(idx) = [];
sd(idx) = [];
x_face = x; m_face = m;sd_face = sd;
subplot(8,4,6);hold on;
myshadedErrorBar(x_body, m_body, sd_body, '-', cols(1, :), 0.4);
myshadedErrorBar(x_face, m_face, sd_face, '-', cols(2, :), 0.4);
line(x_body,m_body,'Color',cols(1, :),'LineWidth',2);
line(x_face,m_face,'Color',cols(2, :),'LineWidth',2);
vertline(0, 'k--');
ylabel('Mean', 'FontSize', 6);title('M2', 'FontSize', 6)
set(gca, 'XLim', [-time_pre_stim, time_per_trial], 'XTick', [-time_pre_stim 0:0.5:time_per_trial], 'YLim', [0 4], 'YTick', 0:4, 'YTickLabel', [], 'FontSize', 5);offsetAxes;
%% sorted variance explained

idx2use = good_units;

% if ~include_eye
% plot the %VE of psth, fig 2C
subplot(8,4,[3,4]);hold on;
[counts, bins] = hist(psth_rsq(idx2use), 0:0.02:1);
histogram(psth_rsq(idx2use), 0:0.02:1, 'FaceColor',  cols(8, :), 'LineStyle', 'none');
set(gca, 'YLim', [0 max(counts)], 'YTick', [0 max(counts)], 'XLim', [bins(min(find(counts == 1))) 1], 'XTick', [bins(min(find(counts == 1))) 1], 'FontSize', 5);vertline(mean(psth_rsq(idx2use)), 'k--');xlabel({'%Variance of SDF', 'explained by model-prediction'}, 'FontSize', 6); ylabel('# units', 'FontSize', 6);offsetAxes;
% end
% plot the %VE for full model, and task model sorted by units
subplot(8,9,[19,20,28,29]);hold on;
[~, idx] = sort(unit_full.all_epochs(idx2use, :), 'descend');
area_mod = unit_task.all_epochs(idx2use, :);
area_full = unit_full.all_epochs(idx2use, :);
plot(area_full(idx), '-', 'Color', cols(3, :), 'LineWidth', 2);
plot(area_mod(idx), '-', 'Color', cols(6, :), 'LineWidth', 1);
axis tight;
set(gca, 'YLim', [-0.05 1], 'YTick', 0:0.5:1, 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
title('All epochs', 'FontSize', 7);ylabel('%VE', 'FontSize', 7);
offsetAxes;

subplot(8,9,[37,38]);hold on;
area_diff = area_full - area_mod;
plot(area_diff(idx), '-', 'Color', [0, 0, 1], 'LineWidth', 0.5);
axis tight;
set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
ylabel('\Delta cv R^{2}', 'FontSize', 7);offsetAxes;

subplot(8,9,39);hold on;
[counts,bins] = hist(area_diff(idx), -0.02:0.01:0.1); %# get counts and bin locations
barh(bins,counts, 'FaceColor', [0, 0, 1], 'LineStyle', 'none');set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'YTickLabel', [], 'XLim', [0 max(counts)], 'XTick', [0 max(counts)], 'FontSize', 6);offsetAxes;

subplot(8,9,[22,23,31,32]);hold on;
[~, idx] = sort(unit_full.controlled_epochs(idx2use, :), 'descend');
area_mod = unit_task.controlled_epochs(idx2use, :);
area_full = unit_full.controlled_epochs(idx2use, :);
plot(area_full(idx), '-', 'Color', cols(3, :), 'LineWidth', 2);
plot(area_mod(idx), '-', 'Color', cols(6, :), 'LineWidth', 1);
axis tight;
set(gca, 'YLim', [-0.05 1], 'YTick', 0:0.5:1, 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
title('Controlled epochs', 'FontSize', 7);
offsetAxes;

subplot(8,9,[40,41]);hold on;
area_diff = area_full - area_mod;
plot(area_diff(idx), '-', 'Color', [0, 0, 1], 'LineWidth', 0.5);
axis tight;
set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
offsetAxes;

subplot(8,9,42);hold on;
[counts,bins] = hist(area_diff(idx), -0.02:0.01:0.1); %# get counts and bin locations
barh(bins,counts, 'FaceColor', [0, 0, 1], 'LineStyle', 'none');set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'YTickLabel', [], 'XLim', [0 max(counts)], 'XTick', [0 max(counts)], 'FontSize', 6);offsetAxes;

subplot(8,9,[25,26,34,35]);hold on;
[~, idx] = sort(unit_full.uncontrolled_epochs(idx2use, :), 'descend');
area_mod = unit_task.uncontrolled_epochs(idx2use, :);
area_full = unit_full.uncontrolled_epochs(idx2use, :);
plot(area_full(idx), '-', 'Color', cols(3, :), 'LineWidth', 2);
plot(area_mod(idx), '-', 'Color', cols(6, :), 'LineWidth', 1);
axis tight;
set(gca, 'YLim', [-0.05 1], 'YTick', 0:0.5:1, 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
title('Uncontrolled epochs', 'FontSize', 7);
offsetAxes;

subplot(8,9,[43,44]);hold on;
area_diff = area_full - area_mod;
plot(area_diff(idx), '-', 'Color', [0, 0, 1], 'LineWidth', 0.5);
axis tight;
set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'XLim', [1 length(idx)], 'XTick', [1 length(idx)], 'FontSize', 6);
offsetAxes;

subplot(8,9,45);hold on;
[counts,bins] = hist(area_diff(idx), -0.02:0.01:0.1); %# get counts and bin locations
barh(bins,counts, 'FaceColor', [0, 0, 1], 'LineStyle', 'none');set(gca, 'YLim', [-0.02 0.1], 'YTick', [-0.02 0:0.05:0.1], 'YTickLabel', [], 'XLim', [0 max(counts)], 'XTick', [0 max(counts)], 'FontSize', 6);offsetAxes;

%% plot unique variance metrics
idx2use = good_units;
ylim_val = [0, .2];
ytick_val = 0:.1:.2;

subplot(4,6,19);hold on;
if ~include_eye
    plot([0.5 3.5], [0 0], 'k--', 'LineWidth', 0.25);
    plot([0.5 3.5], [0.015 0.015], 'k-', 'LineWidth', 0.25);
    y = [all_epochs.task_unique(idx2use), all_epochs.drift_unique(idx2use), all_epochs.move_unique(idx2use)];
    colors2use = cols([6, 8, 3], :);
else
    plot([0.5 5.5], [0 0], 'k--', 'LineWidth', 0.25);
    plot([0.5 5.5], [0.015 0.015], 'k-', 'LineWidth', 0.25);
    y = [all_epochs.task_unique(idx2use), all_epochs.drift_unique(idx2use), all_epochs.move_unique(idx2use), all_epochs.eye_unique(idx2use), all_epochs.pupil_unique(idx2use)];
    colors2use = cols([6, 8, 3, 4, 5], :);
end
for i = 1:size(y,2)
    my_boxPlot(i, y(:, i), colors2use(i, :), 0.5, 0.6);
end
if ~include_eye
    set(gca, 'box', 'off', 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'Task', 'Drift', 'Face+Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
else
    set(gca, 'box', 'off', 'XLim', [0.5, 5.5], 'XTick', 1:5, 'XTickLabel', {'Task', 'Drift', 'Face+Body', 'Eye', 'Pupil'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
end
title(sprintf('All regions, # units = %d', length(idx2use)), 'FontSize', 7);offsetAxes;

subplot(4,6,[20, 21]);hold on;
if ~include_eye
    plot([0.5 9.5], [0 0], 'k--', 'LineWidth', 0.25);
    plot([0.5 9.5], [0.015 0.015], 'k-', 'LineWidth', 0.25);
else
    plot([0.5 15.5], [0 0], 'k--', 'LineWidth', 0.25);
    plot([0.5 15.5], [0.015 0.015], 'k-', 'LineWidth', 0.25);
end
num_units = [];
for i = 1:3
    area_indices = find(sess_units == i);
    idx2use = intersect(good_units, area_indices);
    if ~include_eye
        y = [all_epochs.task_unique(idx2use), all_epochs.drift_unique(idx2use), all_epochs.move_unique(idx2use)];
        colors2use = cols([6, 8, 3], :);
    else
        y = [all_epochs.task_unique(idx2use), all_epochs.drift_unique(idx2use), all_epochs.move_unique(idx2use), all_epochs.eye_unique(idx2use), all_epochs.pupil_unique(idx2use)];
        colors2use = cols([6, 8, 3, 4, 5], :);
    end
    for j = 1:size(y,2)
        if ~include_eye
            my_boxPlot(3*(i-1)+j, y(:, j), colors2use(j, :), 0.5, 0.6);
        else
            my_boxPlot(5*(i-1)+j, y(:, j), colors2use(j, :), 0.5, 0.6);
        end
    end
    num_units = [num_units length(idx2use)];
end
if ~include_eye
    set(gca, 'XLim', [0.5, 9.5], 'XTick', 1:9, 'XTickLabel', {'Task', 'Drift', 'Face+Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
else
    set(gca, 'XLim', [0.5, 15.5], 'XTick', 1:15, 'XTickLabel', {'Task', 'Drift', 'Face+Body', 'Eye', 'Pupil'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
end
title(sprintf('# units: V1 = %d, V2 = %d, V3 = %d', num_units(1), num_units(2), num_units(3)), 'FontSize', 7);offsetAxes;

%% split controlled and uncontrolled retinal inputs
ylim_val = [-0.001, .015];
ytick_val = 0:.005:.015;
subplot(4,6,22);hold on;
idx2use = good_units;
plot([0.3 4.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [controlled_epochs.move_unique(idx2use), controlled_epochs.face_unique(idx2use), controlled_epochs.body_unique(idx2use)];
y2 = [uncontrolled_epochs.move_unique(idx2use), uncontrolled_epochs.face_unique(idx2use), uncontrolled_epochs.body_unique(idx2use)];
colors2use = cols([3, 2, 1], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(controlled_epochs.move_unique(idx2use), uncontrolled_epochs.move_unique(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(controlled_epochs.face_unique(idx2use), uncontrolled_epochs.face_unique(idx2use));
mysigstar(2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(controlled_epochs.body_unique(idx2use), uncontrolled_epochs.body_unique(idx2use));
mysigstar(3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title(sprintf('All regions, # units = %d', length(idx2use)), 'FontSize', 7);offsetAxes;

subplot(4,6,[23, 24]);hold on;
plot([0.3 9.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [controlled_epochs.move_unique(idx2use), controlled_epochs.face_unique(idx2use), controlled_epochs.body_unique(idx2use)];
    y2 = [uncontrolled_epochs.move_unique(idx2use), uncontrolled_epochs.face_unique(idx2use), uncontrolled_epochs.body_unique(idx2use)];
    colors2use = cols([3, 2, 1], :);
    for j = 1:size(y1,2)
        my_boxPlot(3*(i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot(3*(i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    pval = permtest(controlled_epochs.move_unique(idx2use), uncontrolled_epochs.move_unique(idx2use));
    mysigstar(3*(i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(controlled_epochs.face_unique(idx2use), uncontrolled_epochs.face_unique(idx2use));
    mysigstar(3*(i-1)+2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(controlled_epochs.body_unique(idx2use), uncontrolled_epochs.body_unique(idx2use));
    mysigstar(3*(i-1)+3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
end
set(gca, 'XLim', [0.5, 9.5], 'XTick', 1:9, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title(sprintf('# units: V1 = %d, V2 = %d, V3 = %d', num_units(1), num_units(2), num_units(3)), 'FontSize', 7);offsetAxes;

%% polish figures
h=gcf;
set(h, 'PaperType', 'A4');
set(h,'PaperOrientation','portrait');
if ~include_eye
    exportgraphics(gcf,sprintf('%s/fig2_animal_%s.pdf', results_dir, animal),'ContentType','vector');
else
    exportgraphics(gcf,sprintf('%s/fig2_animal_%s_withEye.pdf', results_dir, animal),'ContentType','vector');
end
end

function my_boxPlot(x, y, col, width, patch_alpha)
% plot([x, x], [min(y), max(y)], '_-',  'Color', col, 'LineWidth', 1);
plot([x, x], prctile(y, [100/6, 500/6]), '_-',  'Color', col, 'LineWidth', 1);
patch_x = [x-0.5*width, x+0.5*width, x+0.5*width, x-0.5*width];
plot(patch_x(1:2), [median(y), median(y)], '-',  'Color', col, 'LineWidth', 1);
iqr_y = prctile(y, [25, 75]);
patch_y = [iqr_y(1) iqr_y(1) iqr_y(2) iqr_y(2)];
patch(patch_x, patch_y, col,'FaceColor', col, 'EdgeColor', col, 'FaceAlpha', patch_alpha);
end