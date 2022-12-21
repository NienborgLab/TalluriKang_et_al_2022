function [mean_rates, mean_quartile_rates] = unit_activity(Vc, D)

[~, stim_idx] = intersect(D.timepts2use, D.stim_timepts);
mean_rates = sum(Vc(:, stim_idx), 2)/(D.dt*(length(find(diff(D.stim_timepts) == 1))-1));

% compute mean stim rates for each 25th percentile of session
quartile_breaks = prctile(1:D.timepts2use(end), [0 25 50 75 100]);
mean_quartile_rates = nan(size(Vc, 1), 4);

for i = 1:4
    [~, stim_idx, temp_idx] = intersect(D.timepts2use(find(D.timepts2use >= quartile_breaks(i) & D.timepts2use < quartile_breaks(i+1))), D.stim_timepts);
    mean_quartile_rates(:, i) = sum(Vc(:, stim_idx), 2)/(D.dt*(length(find(diff(D.stim_timepts(temp_idx)) == 1))-1));
end
end