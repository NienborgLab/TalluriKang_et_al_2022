function [nonStimEvents, nonStimEventType, nonStimLabels] = getNonTaskVariabels(D)
% Bharath Talluri, Dec 2022

% Get NonTask Variables will set up the any non-stimulus events for
% processing by the design matrix

% Inputs:
%   D <struct>: The preprocessed data structure (see code/preprocess/M2/preprocess_dataset)

% Outputs;
%   nonStimEvents    <double>: NT x NCovariates, the events themselves
%   nonStimEventType <double>: which event type they are (1,2, or 3)
%   nonStimLabels      <cell>: cell array of covariate labels

% animal M2 was presented with two targets after stimulus presentation.
% This was followed by the monkey making a saccade towards the chosen
% alternative, followed by a reward if the choice was correct.
% So we have six non-stimulus events in a trial: time, target locations, 
% motor response, choice and reward. Time regressor models fluctuations 
% as a function of time for the duration of the trial starting at the onset
% of stimulus presentation and target and reward regressors model 
% fluctuations due to target and reward presentation respectively.These are
% peri-event type variables. Choice and motor response regressors also
% model preparatory activity for choice and saccades that most likely
% precede the event. 

nonStimLabels = {'time', 'targ_left', 'targ_right', 'motor_response', 'choice', 'reward'}; %some non-stimulus variables, could include microsaccades in the future..
nonStimEventType = [1, 2, 2, 3, 3, 2]; %different type of events. these are all peri-event variables.
nonStimEvents(:,1) = D.stim_start; % event type 1; label 'time'
nonStimEvents(:,2) = D.target_left;
nonStimEvents(:,3) = D.target_right;
nonStimEvents(:,4) = D.resp;
nonStimEvents(:,5) = D.choice;
nonStimEvents(:,6) = D.reward_start;