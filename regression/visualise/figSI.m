function figSI(varargin)

ip = inputParser();
ip.addParameter('include_eye', 0)
ip.parse(varargin{:});

include_eye = ip.Results.include_eye;

% define directories
results_dir = '../data/results';

no_eye = load(sprintf('%s/harvested_vars.mat', results_dir));
with_eye = load(sprintf('%s/harvested_vars_withEye.mat', results_dir));

figure(1);clf;
colormap linspecer;
cols = linspecer(9,'qualitative');
cols(1, :) = [0 0 1];
cols(2, :) = [0 0.65 1];

load(sprintf('%s/selected_units.mat', results_dir));

%% compare with and without eye
ylim_val = [-0.0005, .002];
ytick_val = 0:.001:.002;

subplot(4,6,1);hold on;
idx2use = good_units;
plot([0.3 4.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [no_eye.controlled_epochs.move_unique(idx2use), no_eye.controlled_epochs.face_unique(idx2use), no_eye.controlled_epochs.body_unique(idx2use)];
y2 = [with_eye.controlled_epochs.move_unique(idx2use), with_eye.controlled_epochs.face_unique(idx2use), with_eye.controlled_epochs.body_unique(idx2use)];
colors2use = cols([3, 2, 1], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(no_eye.controlled_epochs.move_unique(idx2use), with_eye.controlled_epochs.move_unique(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(no_eye.controlled_epochs.face_unique(idx2use), with_eye.controlled_epochs.face_unique(idx2use));
mysigstar(2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(no_eye.controlled_epochs.body_unique(idx2use), with_eye.controlled_epochs.body_unique(idx2use));
mysigstar(3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
offsetAxes;

subplot(4,6,[2, 3]);hold on;
plot([0.3 9.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(no_eye.sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [no_eye.controlled_epochs.move_unique(idx2use), no_eye.controlled_epochs.face_unique(idx2use), no_eye.controlled_epochs.body_unique(idx2use)];
    y2 = [with_eye.controlled_epochs.move_unique(idx2use), with_eye.controlled_epochs.face_unique(idx2use), with_eye.controlled_epochs.body_unique(idx2use)];
    colors2use = cols([3, 2, 1], :);
    for j = 1:size(y1,2)
        my_boxPlot(3*(i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot(3*(i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    pval = permtest(no_eye.controlled_epochs.move_unique(idx2use), with_eye.controlled_epochs.move_unique(idx2use));
    mysigstar(3*(i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(no_eye.controlled_epochs.face_unique(idx2use), with_eye.controlled_epochs.face_unique(idx2use));
    mysigstar(3*(i-1)+2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(no_eye.controlled_epochs.body_unique(idx2use), with_eye.controlled_epochs.body_unique(idx2use));
    mysigstar(3*(i-1)+3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
end
set(gca, 'XLim', [0.5, 9.5], 'XTick', 1:9, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('Controlled input epochs', 'FontSize', 7);offsetAxes;

ylim_val = [-0.0015, .015];
ytick_val = 0:.005:.015;

subplot(4,6,4);hold on;
idx2use = good_units;
plot([0.3 4.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [no_eye.uncontrolled_epochs.move_unique(idx2use), no_eye.uncontrolled_epochs.face_unique(idx2use), no_eye.uncontrolled_epochs.body_unique(idx2use)];
y2 = [with_eye.uncontrolled_epochs.move_unique(idx2use), with_eye.uncontrolled_epochs.face_unique(idx2use), with_eye.uncontrolled_epochs.body_unique(idx2use)];
colors2use = cols([3, 2, 1], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(no_eye.uncontrolled_epochs.move_unique(idx2use), with_eye.uncontrolled_epochs.move_unique(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(no_eye.uncontrolled_epochs.face_unique(idx2use), with_eye.uncontrolled_epochs.face_unique(idx2use));
mysigstar(2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(no_eye.uncontrolled_epochs.body_unique(idx2use), with_eye.uncontrolled_epochs.body_unique(idx2use));
mysigstar(3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
offsetAxes;

subplot(4,6,[5, 6]);hold on;
plot([0.3 9.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(no_eye.sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [no_eye.uncontrolled_epochs.move_unique(idx2use), no_eye.uncontrolled_epochs.face_unique(idx2use), no_eye.uncontrolled_epochs.body_unique(idx2use)];
    y2 = [with_eye.uncontrolled_epochs.move_unique(idx2use), with_eye.uncontrolled_epochs.face_unique(idx2use), with_eye.uncontrolled_epochs.body_unique(idx2use)];
    colors2use = cols([3, 2, 1], :);
    for j = 1:size(y1,2)
        my_boxPlot(3*(i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot(3*(i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    pval = permtest(no_eye.uncontrolled_epochs.move_unique(idx2use), with_eye.uncontrolled_epochs.move_unique(idx2use));
    mysigstar(3*(i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(no_eye.uncontrolled_epochs.face_unique(idx2use), with_eye.uncontrolled_epochs.face_unique(idx2use));
    mysigstar(3*(i-1)+2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(no_eye.uncontrolled_epochs.body_unique(idx2use), with_eye.uncontrolled_epochs.body_unique(idx2use));
    mysigstar(3*(i-1)+3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
end
set(gca, 'XLim', [0.5, 9.5], 'XTick', 1:9, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('Uncontrolled input epochs', 'FontSize', 7);offsetAxes;

clear no_eye with_eye good_units
%% split inferred and uninferred retinal inputs during uncontrolled retinal inout epochs

if ~include_eye
    load(sprintf('%s/harvested_vars.mat', results_dir));
    ylim_val = [-0.005, .016];
    ytick_val = -0.005:.005:.015;
else
    load(sprintf('%s/harvested_vars_withEye.mat', results_dir));
    ylim_val = [-0.009, .006];
    ytick_val = -0.006:.006:.006;
end
load(sprintf('%s/selected_units.mat', results_dir));
good_units = intersect(good_units, find(~isnan(uncontrolled_epochs.infer.move_unique) & ~isnan(uncontrolled_epochs.uninfer.move_unique)));

subplot(4,6,7);hold on;
idx2use = good_units;
plot([0.3 4.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [uncontrolled_epochs.infer.move_unique(idx2use), uncontrolled_epochs.infer.face_unique(idx2use), uncontrolled_epochs.infer.body_unique(idx2use)];
y2 = [uncontrolled_epochs.uninfer.move_unique(idx2use), uncontrolled_epochs.uninfer.face_unique(idx2use), uncontrolled_epochs.uninfer.body_unique(idx2use)];
colors2use = cols([3, 2, 1], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(uncontrolled_epochs.infer.move_unique(idx2use), uncontrolled_epochs.uninfer.move_unique(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(uncontrolled_epochs.infer.face_unique(idx2use), uncontrolled_epochs.uninfer.face_unique(idx2use));
mysigstar(2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
pval = permtest(uncontrolled_epochs.infer.body_unique(idx2use), uncontrolled_epochs.uninfer.body_unique(idx2use));
mysigstar(3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title(sprintf('All regions, # units = %d', length(idx2use)), 'FontSize', 7);offsetAxes;

num_units = [];
subplot(4,6,[8, 9]);hold on;
plot([0.3 9.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [uncontrolled_epochs.infer.move_unique(idx2use), uncontrolled_epochs.infer.face_unique(idx2use), uncontrolled_epochs.infer.body_unique(idx2use)];
    y2 = [uncontrolled_epochs.uninfer.move_unique(idx2use), uncontrolled_epochs.uninfer.face_unique(idx2use), uncontrolled_epochs.uninfer.body_unique(idx2use)];
    colors2use = cols([3, 2, 1], :);
    for j = 1:size(y1,2)
        my_boxPlot(3*(i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot(3*(i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    pval = permtest(uncontrolled_epochs.infer.move_unique(idx2use), uncontrolled_epochs.uninfer.move_unique(idx2use));
    mysigstar(3*(i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(uncontrolled_epochs.infer.face_unique(idx2use), uncontrolled_epochs.uninfer.face_unique(idx2use));
    mysigstar(3*(i-1)+2, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    pval = permtest(uncontrolled_epochs.infer.body_unique(idx2use), uncontrolled_epochs.uninfer.body_unique(idx2use));
    mysigstar(3*(i-1)+3, max(ylim_val)-0.1*max(ylim_val), pval, 0);
    num_units = [num_units length(idx2use)];
end
set(gca, 'XLim', [0.5, 9.5], 'XTick', 1:9, 'XTickLabel', {'Face+Body', 'Face', 'Body'}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title(sprintf('# units: V1 = %d, V2 = %d, V3 = %d', num_units(1), num_units(2), num_units(3)), 'FontSize', 7);offsetAxes;

%% compare high and low movement trials
load(sprintf('%s/selected_units.mat', results_dir));
if ~include_eye
    ylim_val = [-0.002, .002];
    ytick_val = -0.002:.002:.002;
else
    ylim_val = [-0.002, .002];
    ytick_val = -0.002:.002:.002;
end

subplot(4,6,13);hold on;
idx2use = good_units;
y1 = [face_var_trl_split.mean.low face_var_trl_split.var.low];
y2 = [face_var_trl_split.mean.high face_var_trl_split.var.high];
colors2use = cols([3, 3], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(face_var_trl_split.mean.high, face_var_trl_split.mean.low);
mysigstar(1, 1, pval, 0);
pval = permtest(face_var_trl_split.var.high, face_var_trl_split.var.low);
mysigstar(2, 1.4, pval, 0);
set(gca, 'XLim', [0.5, 2.5], 'XTick', 1:2, 'XTickLabel', {'Mean', 'Variance'}, 'XTickLabelRotation', -45, 'YLim', [0 1.5], 'YTick', 0:0.5:1.5, 'FontSize', 6);
title('Face Camera', 'FontSize', 7);offsetAxes;

subplot(4,6,14);hold on;
idx2use = good_units;
plot([0.3 1.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [face_var_trl_split.rsq.low(idx2use)];
y2 = [face_var_trl_split.rsq.high(idx2use)];
colors2use = cols(3, :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
plot([0.3 1.7], [mean(controlled_epochs.move_unique(idx2use)) mean(controlled_epochs.move_unique(idx2use))], 'k', 'LineWidth', 1);
pval = permtest(face_var_trl_split.rsq.high(idx2use), face_var_trl_split.rsq.low(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 1.5], 'XTick', [], 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('High vs Low', 'FontSize', 7);%offsetAxes;

subplot(4,6,[15, 16]);hold on;
plot([0.3 3.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [face_var_trl_split.rsq.low(idx2use)];
    y2 = [face_var_trl_split.rsq.high(idx2use)];
    for j = 1:size(y1,2)
        my_boxPlot((i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot((i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    plot([(i-1)+0.6 (i-1)+1.4], [mean(controlled_epochs.move_unique(idx2use)) mean(controlled_epochs.move_unique(idx2use))], 'k', 'LineWidth', 1);
    pval = permtest(face_var_trl_split.rsq.high(idx2use), face_var_trl_split.rsq.low(idx2use));
    mysigstar((i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
end
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'', ''}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('Face Movement trials', 'FontSize', 7);offsetAxes;

subplot(4,6,19);hold on;
idx2use = good_units;
y1 = [body_var_trl_split.mean.low body_var_trl_split.var.low];
y2 = [body_var_trl_split.mean.high body_var_trl_split.var.high];
colors2use = cols([3, 3], :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
pval = permtest(body_var_trl_split.mean.high, body_var_trl_split.mean.low);
mysigstar(1, 1, pval, 0);
pval = permtest(body_var_trl_split.var.high, body_var_trl_split.var.low);
mysigstar(2, 1.4, pval, 0);
set(gca, 'XLim', [0.5, 2.5], 'XTick', 1:2, 'XTickLabel', {'Mean', 'Variance'}, 'XTickLabelRotation', -45, 'YLim', [0 2.5], 'YTick', 0:0.5:2.5, 'FontSize', 6);
title('Body Camera', 'FontSize', 7); offsetAxes;

subplot(4,6,20);hold on;
idx2use = good_units;
plot([0.3 1.7], [0 0], 'k--', 'LineWidth', 0.25);
y1 = [body_var_trl_split.rsq.low(idx2use)];
y2 = [body_var_trl_split.rsq.high(idx2use)];
colors2use = cols(3, :);
for i = 1:size(y1,2)
    my_boxPlot(i-0.2, y1(:, i), colors2use(i, :), 0.4, 0.6);
    my_boxPlot(i+0.2, y2(:, i), colors2use(i, :), 0.4, 0);
end
plot([0.3 1.7], [mean(controlled_epochs.move_unique(idx2use)) mean(controlled_epochs.move_unique(idx2use))], 'k', 'LineWidth', 1);
pval = permtest(face_var_trl_split.rsq.high(idx2use), face_var_trl_split.rsq.low(idx2use));
mysigstar(1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
set(gca, 'XLim', [0.5, 1.5], 'XTick', [], 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('High vs Low', 'FontSize', 7);%offsetAxes;


subplot(4,6,[21, 22]);hold on;
plot([0.3 3.7], [0 0], 'k-', 'LineWidth', 0.25);
for i = 1:3
    area_indices = find(sess_units == i);
    idx2use = intersect(good_units, area_indices);
    y1 = [body_var_trl_split.rsq.low(idx2use)];
    y2 = [body_var_trl_split.rsq.high(idx2use)];
    for j = 1:size(y1,2)
        my_boxPlot((i-1)+j-0.2, y1(:, j), colors2use(j, :), 0.4, 0.6);
        my_boxPlot((i-1)+j+0.2, y2(:, j), colors2use(j, :), 0.4, 0);
    end
    plot([(i-1)+0.6 (i-1)+1.4], [mean(controlled_epochs.move_unique(idx2use)) mean(controlled_epochs.move_unique(idx2use))], 'k', 'LineWidth', 1);
    pval = permtest(face_var_trl_split.rsq.high(idx2use), face_var_trl_split.rsq.low(idx2use));
    mysigstar((i-1)+1, max(ylim_val)-0.1*max(ylim_val), pval, 0);
end
set(gca, 'XLim', [0.5, 3.5], 'XTick', 1:3, 'XTickLabel', {'', ''}, 'XTickLabelRotation', -45, 'YLim', ylim_val, 'YTick', ytick_val, 'YTickLabel', 100*ytick_val, 'FontSize', 6);
title('Body Movement trials', 'FontSize', 7);offsetAxes;

%% polish figures
h=gcf;
set(h, 'PaperType', 'A4');
set(h,'PaperOrientation','portrait');
if ~include_eye
    exportgraphics(gcf,sprintf('%s/figSI.pdf', results_dir),'ContentType','vector');
else
    exportgraphics(gcf,sprintf('%s/figSI_withEye.pdf', results_dir),'ContentType','vector');
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