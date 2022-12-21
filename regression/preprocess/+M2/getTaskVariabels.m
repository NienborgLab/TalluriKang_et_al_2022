function [taskEvents, taskEventType, taskLabels] = getTaskVariabels(D, varargin)
% Written by J. Yates, 2021
% Modified for sharing by Bharath Talluri, Dec 2022

% Get Task Variables will set up the task events for processing by the
% design matrix

% Inputs:
%   D <struct>: The data struct (see io.preprocess_data)

% Optional input (as argument pairs):
%   'attention_mode': ['none', 'full'(default), 'signed']
%   'include_contrast': true or false (default: true)

% Outputs;
%   taskEvents    <double>: NT x NCovariates, the events themselves
%   taskEventType <double>: which event type they are (1,2, or 3)
%   taskLabels      <cell>: cell array of covariate labels

ip = inputParser();
ip.addParameter('attention_mode', 'full')
ip.addParameter('include_contrast', true)
ip.parse(varargin{:});

% Unpack the disparity stimulus so each value is represented
hdx_left_labels = arrayfun(@(x) sprintf('hdx_left_%02.2f', x), D.hdx_levels, 'uni', 0);
hdx_right_labels = arrayfun(@(x) sprintf('hdx_right_%02.2f', x), D.hdx_levels, 'uni', 0);

nhdx = numel(D.hdx_levels);

clear taskEvents
taskLabels = [];
taskEventType = [];
ctr = 0; % keep track of index for covariate

if strcmp(ip.Results.attention_mode, 'full')

    % we can include the full potential of attention to change tuning by
    % simply separating the covariates for the attend IN and OUT conditions
    % This doubles the number of stimulus parameters. These are all 
    % peri-event types

    % covariate 1: hdx_left & attend IN
    taskLabels = [taskLabels cellfun(@(x) ['IN_' x], hdx_left_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.hdx_left) .* D.attend) < 0;
    ctr = ctr + nhdx;
    
    % covariate 2: hdx_left & attend OUT
    taskLabels = [taskLabels cellfun(@(x) ['OUT_' x], hdx_left_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.hdx_left) .* D.attend) > 0;
    ctr = ctr + nhdx;

    % covariate 3: hdx_right & attend IN
    taskLabels = [taskLabels cellfun(@(x) ['IN_' x], hdx_right_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.hdx_right) .* D.attend) > 0;
    ctr = ctr + nhdx;

    % covariate 4: hdx_right & attend OUT
    taskLabels = [taskLabels cellfun(@(x) ['OUT_' x], hdx_right_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.hdx_right) .* D.attend) < 0;
    ctr = ctr + nhdx;

else

    % here, no matter how we handle attention, the stimulus covariates
    % are untouched by it

    % covariate 1: hdx_left
    taskLabels = [taskLabels hdx_left_labels]; %#ok<*UNRCH>
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = D.hdx_left;
    ctr = ctr + nhdx;

    % covariate 2: hdx_right
    taskLabels = [taskLabels hdx_right_labels];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = D.hdx_right;
    ctr = ctr + nhdx;

    if strcmp(ip.Results.attention_mode, 'signed')
        % here, we treat attention as another additive variable that starts
        % when contrast is on
        % covariate: Attention
        % (eventType 1, means it effects for the rest of the trial)
        stim_onset = [0; diff(D.contrast_left)==1];
        taskEvents(:,ctr+1) = stim_onset .* D.attend;
        taskLabels = [taskLabels 'Attention'];
        taskEventType = [taskEventType 1];
        ctr = ctr + 1;
    end

end

if ip.Results.include_contrast
    % --- "Contrast" covariates add the stimulus onset
    % covariate: contrast left
    taskEvents(:,ctr+1) = [0; diff(D.contrast_left)==1];
    taskLabels = [taskLabels 'contrast_left'];
    taskEventType = [taskEventType 1];
    ctr = ctr + 1;

    % covariate: contrast right
    taskEvents(:,ctr+1) = [0; diff(D.contrast_right)==1];
    taskLabels = [taskLabels 'contrast_right'];
    taskEventType = [taskEventType 1];
    ctr = ctr + 1;
end