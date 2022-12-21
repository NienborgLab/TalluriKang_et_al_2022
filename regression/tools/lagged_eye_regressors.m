function lagged_regressors = lagged_eye_regressors(regressor, opts, blocks)

% Bharath Talluri, January 2022
% convolve video regressors with tent basis functions
% adapted from code by Stefano Panzeri's group in Runyan et al.

% regressor is time x 1
% output is time x nbasis
% define the basis functions first
dt = opts.dt;
twin = opts.eye_twin;

% define basis functions first
basis_centers = -twin:twin/2:twin/2;nBasis = length(basis_centers);
NT=-twin:dt:twin;step = ceil((length(NT)-1)/nBasis);basis_func = max(1 - abs( (NT)'-basis_centers)/(dt*step), 0);
bases = basis_func(:, 2:end)';

% convolution potentially leads to edge artefacts: https://github.com/distillpub/post--deconv-checkerboard/issues/5
% so we can do symmetric padding: copying the timeseries at the ends
% the length of symmetric padding depends on the length of the basis
% functions
padding = size(bases,2);
trialCnt = size(blocks,1); %nr of trials
lagged_regressors = nan(size(regressor, 1), size(bases, 1));
% convolve regressor with individual basis functions
for b = 1:size(bases, 1)
    for t = 1:trialCnt
        % do not max-normalise the convolved regressor because some
        % parametric regressors might take nan values. Also,
        % max-normalisation is not required here since we standardise
        % individual regressors before fitting.
        this_trial_regressor = regressor(blocks(t, 1):blocks(t, 2));
        % too many nans in trial regresor will cause issue with convolution
        this_trial_regressor = fillmissing(this_trial_regressor, 'movmedian', 10);
        padded_regressor = [this_trial_regressor(padding:-1:1);this_trial_regressor;this_trial_regressor(end:-1:end-padding+1)]';
        convolved_regressor = nanconv(padded_regressor, bases(b, :), 'noedge', '1d');
        lagged_regressors(blocks(t, 1):blocks(t, 2), b) = convolved_regressor(1+padding:end-padding);
        clear this_trial_regressor padded_regressor convolved_regressor
    end
end
end