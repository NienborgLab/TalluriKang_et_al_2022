function [taskEvents, taskEventType, taskLabels] = getTaskVariabels(D, varargin)
% Written by J. Yates, 2021
% Modified for sharing by Bharath Talluri, Dec 2022

% Get Task Variables will set up the task events for processing by the
% design matrix

% Inputs:
%   D <struct>: The data struct (see io.preprocess_data)

% Optional input (as argument pairs):
%   'position_mode': ['none', 'full' (default)]

% Outputs;
%   taskEvents    <double>: NT x NCovariates, the events themselves
%   taskEventType <double>: which event type they are (1,2, or 3)
%   taskLabels      <cell>: cell array of covariate labels

ip = inputParser();
ip.addParameter('position_mode', 'full')
ip.parse(varargin{:});

% Unpack the stimulus so each contrast value (but 0) is represented
cont_labels = arrayfun(@(x) sprintf('hdx_%4.4f', x), D.cont_levels(D.cont_levels ~= 1001), 'uni', 0);
nhdx = numel(cont_labels);

clear taskEvents
taskLabels = [];
taskEventType = [];
ctr = 0; % keep track of index for covariate

if strcmp(ip.Results.position_mode, 'full')
    % in this mode, we will have separate regressors for each position, and
    % each contrast value to account for any variability as a function of
    % the sample position in the trial. All stimuli are peri-event type
    % variables

    % covariate 1: sample position 1
    taskLabels = [taskLabels cellfun(@(x) [x '_smp1'], cont_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.cont_seq(:, 1:nhdx)) .* double(D.smp1)) == 1;
    ctr = ctr + nhdx;

    % covariate 2: sample position 2
    taskLabels = [taskLabels cellfun(@(x) [x '_smp2'], cont_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.cont_seq(:, 1:nhdx)) .* double(D.smp2)) == 1;
    ctr = ctr + nhdx;

    % covariate 3: sample position 3
    taskLabels = [taskLabels cellfun(@(x) [x '_smp3'], cont_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.cont_seq(:, 1:nhdx)) .* double(D.smp3)) == 1;
    ctr = ctr + nhdx;

    % covariate 4: sample position 4
    taskLabels = [taskLabels cellfun(@(x) [x '_smp4'], cont_labels, 'uni', 0)];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = (double(D.cont_seq(:, 1:nhdx)) .* double(D.smp4)) == 1;
    ctr = ctr + nhdx;

elseif strcmp(ip.Results.position_mode, 'none')
    % in this mode, we ignore any variability due to position and only
    % consider variability due to changes in stimulus contrast

    taskLabels = [taskLabels cont_labels];
    taskEventType = [taskEventType 2*ones(1, nhdx)];
    taskEvents(:,ctr+(1:nhdx)) = double(D.cont_seq(:, 1:end-1));
    ctr = ctr + nhdx;
end
end