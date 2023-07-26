function [mean_stim_rates, blockwise_mean_stim_rates] = stim_rates_ma(Vc, D)
% stim_rates = nan(size(Vc, 1), size(D.blocks_clean, 1));
% k = 1;dt = 1/D.frameRate;
% for t = 1:size(D.blocks_clean, 1)
%     stim_start = D.blocks_clean(t, 1) + 0.3/dt;
%     if stim_start + 120 > D.blocks_clean(t, 2)
%         keyboard
%     end
%     stim_rates(:, k) = sum(Vc(:, stim_start: stim_start+120), 2)/2;
%     k = k + 1;
% end
% mean_stim_rates = mean(stim_rates, 2);

[~, stim_idx] = intersect(D.timepts2use, D.stim_timepts);
mean_stim_rates = sum(Vc(:, stim_idx), 2)/(D.dt*(length(find(diff(D.stim_timepts) == 1))-1));

% compute mean stim rates for each 25th percentile of session
block_breaks = prctile(1:D.timepts2use(end), [0 25 50 75 100]);
blockwise_mean_stim_rates = nan(size(Vc, 1), 4);

for i = 1:4
    [~, stim_idx, temp_idx] = intersect(D.timepts2use(find(D.timepts2use >= block_breaks(i) & D.timepts2use < block_breaks(i+1))), D.stim_timepts);
    blockwise_mean_stim_rates(:, i) = sum(Vc(:, stim_idx), 2)/(D.dt*(length(find(diff(D.stim_timepts(temp_idx)) == 1))-1));
end
end