function movement_regressors = lagged_video_regressors(regressor, opts)

% Bharath Talluri, January 2022
% convolve video regressors with tent basis functions
% adapted from code by Stefano Panzeri's group in Runyan et al.

% regressor is time x num_regressors
% output is n_basis x time x num_regressors
% define the basis functions first
dt = opts.dt;
twin = opts.twin;
frames = opts.framesPerTrial;


% define basis functions first
basis_centers = -twin:twin/2:twin/2;nBasis = length(basis_centers);
NT=-twin:dt:twin;step = ceil((length(NT)-1)/nBasis);basis_func = max(1 - abs( (NT)'-basis_centers)/(dt*step), 0);
bases = basis_func(:, 2:end)';

% convolution potentially leads to edge artefacts: https://github.com/distillpub/post--deconv-checkerboard/issues/5
% so we can do symmetric padding: copying the timeseries at the ends
% the length of symmetric padding depends on the length of the basis
% functions
padding = size(bases,2);
regressor = reshape(regressor, frames, [], size(regressor, 2)); %reshape to trials
trialCnt = size(regressor,2); %nr of trials

for iReg = 1:size(regressor, 3)
    this_reg = regressor(:, :, iReg);
    this_reg = this_reg';
    padded_regressor = cat(2,squeeze(this_reg(:,padding:-1:1)),squeeze(this_reg(:,:)),squeeze(this_reg(:,end:-1:end-padding+1)));
    % ntrials = size(this_reg,1);
    % convolve regressor with individual basis functions
    for b = 1:size(bases, 1)
        for t = 1:trialCnt
            % do not max-normalise the convolved regressor because some
            % parametric regressors might take nan values. Also,
            % max-normalisation is not required here since we standardise
            % individual regressors before fitting.
            this_basis_reg_conv(:, t) = conv(padded_regressor(t, :), bases(b, :), 'same');
        end
        movement_regressors(b, :, iReg) = reshape(this_basis_reg_conv(1+padding:end-padding,:), [], 1);
        clear this_basis_reg_conv
    end
end
end