function [nonStimEvents, nonStimEventType, nonStimLabels] = getNonTaskVariabels(D)
% Bharath Talluri, Dec 2022

% Get NonTask Variables will set up the any non-stimulus events for
% processing by the design matrix

% Inputs:
%   D <struct>: The preprocessed data structure (see code/preprocess/M1/preprocess_dataset)

% Outputs;
%   nonStimEvents    <double>: NT x NCovariates, the events themselves
%   nonStimEventType <double>: which event type they are (1,2, or 3)
%   nonStimLabels      <cell>: cell array of covariate labels

% animal M1 got a reward after every trial where it did not break fixation
% during the stimulus presentation window. So we have two non-stimulus
% events in a trial: time, and reward. time regressor models fluctuations 
% as a function of time for the duration of the trial starting at the onset
% of stimulus presentation and reward regressor models fluctuations due to
% reward presentation. Both are peri-event type variables

nonStimLabels = {'time', 'reward'};
nonStimEventType = [1, 2];
nonStimEvents(:,1) = D.stim_start;
nonStimEvents(:,2) = D.reward_start;