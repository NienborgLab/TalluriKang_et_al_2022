function [D, D_eye] = preprocess_dataset(datadir, exname, varargin)
% first written by J. Yates, 2021
% modified by Bharath Talluri, Dec 2021
% Modified for sharing by Bharath Talluri, Dec 2022

% load and preprocess raw data

% Preprocess dataset will take all of the files in an experiment folder,
% time-align them for each trial and then bin and make it time-continuous

% Inputs:
%   datadir <char>: path to experiment director
%   exname  <char>: experiment name (e.g., M1_0008)

% Optional Inputs (as argument pairs)
%   'dt' <double>: bin size (in seconds)
%   'time_pre_stim'  <double>: time (sec) to include before stimulus onset
%   'time_per_trial' <double>: time (sec), duration of each trial

% Outputs:
% D <struct>: time-continuous structure
%   hdx_levels    <double>: the unique disparity levels of the hdx_sequences
%   hdx_right      <uint8>: NT x NLevels, hdx sequence on the right of the
%                           screen, as a "one hot" sequence
%   hdx_left       <uint8>: NT x NLevels, hdx sequence on the left of screen
%   contrast_left  <uint8>: NT x 1, whether (left) stimulus was on
%   contrast_right <uint8>: NT x 1, whether (right) stimulus was on
%   motenergy      <uint8>: NT x 1, motion energy at each time bin
%   attend        <double>: NT x 1, attention condition (-1 or 1, left or right)
%   robs           <uint8>: NT x NCells, binned spike counts
%   blocks        <double>: NTrials x 2, trial start / stop
%   resp            <int8>: NT x 1, signed response at time of response
%   correct       <double>: NT x 1, whether animal was correct (at time of
%                           response)
%   faceSVD        <int16>: NT x 1000, SVD components of face motion energy
%   bodySVD        <int16>: NT x 1000, SVD components of body motion energy
%   units         <strcut>: struct array of unit information
%   frameRate     <double>: time resolution

ip = inputParser();
ip.addParameter('dt', 1/60)
ip.addParameter('time_pre_stim', .3)
ip.addParameter('time_per_trial', 5)
ip.parse(varargin{:});

% load data
fl = dir(fullfile(datadir, [exname '*']));
fpath = fullfile(fl.folder, fl.name);

flist = {dir(fullfile(fpath, '*.mat')).name};
for f = 1:numel(flist)
    if ~isempty(flist{f})
        if strcmp(flist{f}(1:2), '._')
            continue
        end
        if strcmp(flist{f}, 'faceSVD.mat')
            fprintf('Loading [%s]\n', flist{f})
            tmp = load(fullfile(fpath, flist{f}));
            faceSVD = tmp.vSVD;
            clear tmp
        elseif strcmp(flist{f}, 'bodySVD.mat')
            fprintf('Loading [%s]\n', flist{f})
            tmp = load(fullfile(fpath, flist{f}));
            bodySVD = tmp.vSVD;
            clear tmp
        else
            fprintf('Loading [%s]\n', flist{f})
            load(fullfile(fpath, flist{f}));
        end
    end
end

% get spike information
NC = numel(stKs); %#ok<*USENS>
st = [];
clu = [];
for cc = 1:NC
    st_ = stKs{cc}(:);
    st = [st; st_]; %#ok<*AGROW>
    clu = [clu; cc*ones(numel(st_),1)];
end
[st, ind] = sort(st);
clu = clu(ind);

% analyse only successful trials: when the animals fixated through all 4
% stimuli and obtained a reward at the end. Each preprocessed trial thus
% contains 4 samples
full_trials = find([ex.Trials.RewardSize]~=0);
full_trials = full_trials(full_trials >3);
% one session did not have video files spanning the whole session so we remove
% the last few samples where we did not have video data
if strcmp(exname, 'M1_0032')
    full_trials = setdiff(full_trials, length(ex.Trials)-38:length(ex.Trials));
end

nTrials = numel(full_trials);

% save information relevant to the eye regressors
dpp = atan(ex.setup.monitorWidth/2/ex.setup.viewingDistance)*180/pi/(ex.setup.screenRect(3)/2);
if isempty(ex.eyeCal.Delta(1).TrialNo)
    ex.eyeCal.Delta(1).TrialNo = 1;
end
RX0 = ex.eyeCal.Delta(1).RX0;
RY0 = ex.eyeCal.Delta(1).RY0;
LX0 = ex.eyeCal.Delta(1).LX0;
LY0 = ex.eyeCal.Delta(1).LY0;
screen_height = ex.setup.screenRect(4)/2;
screen_width = ex.setup.screenRect(3)/2;

% S is a struct array of length N Trials, and accomplishes time-aligning
% everything within each trial

S = repmat(struct('cont', [], ...
    'cont_time', [], ...
    'cont_seq', [], ...
    'reward_size', [], ...
    'reward_offset', [], ...
    'reward_onset', [], ...
    'fix_start', [], ...
    'fix_lag', [], ...
    'isi_lag', [], ...
    'nframes', [], ...
    'stim_start', [], ...
    'break_fix', [], ...
    'smp1', [], ...
    'smp2', [], ...
    'smp3', [], ...
    'smp4', [], ...
    'eyevals', []), nTrials, 1);

fprintf('Get stimulus information \n');

for iTrial = 1:nTrials
    fprintf('Trial %d/%d\n', iTrial, nTrials)

    trl_num_idx = full_trials(iTrial);
    trl_smp_idx = trl_num_idx-3:trl_num_idx;
    % deal with some session-specific exceptions
    if strcmp(exname, 'M1_0033')
        if ex.Trials(trl_num_idx).x0 ~= -0.25
            continue
        end
    end
    if strcmp(exname, 'M1_0008')
        if ex.Trials(trl_num_idx).x0 ~= -7.0
            continue
        end
    end
    if strcmp(exname, 'M1_0012')
        if ex.Trials(trl_num_idx).x0 ~= -3.5
            continue
        end
    end
    if strcmp(exname, 'M1_0042')
        if ex.Trials(trl_num_idx).y0 ~= 3.6
            continue
        end
    end

    % double-check to make sure we are not including any unrewarded trials
    if any([ex.Trials(trl_smp_idx).Reward] == 0)
        continue
    end

    % start time-alignment within trials
    S(iTrial).sample_onset = [];
    % clock alignment: there are different clocks in the experimental
    % session: the clock from stimulus computer, the clock on the neural recording
    % computer, and the clock on the eye-tracker computer. We need to make
    % sure these clocks are aligned before we can align spikes, stimulus
    % events, and eye/pupil timeseries
    % align the clocks and save the sample onset times in the data
    % structure
    for smp_idx = trl_smp_idx
        ephys_start = ex.Trials(smp_idx).nevTrialStart;
        if strcmp(exname, 'M1_0037')
            ptb_start = ex.Trials(smp_idx).TrialStart;
        else
            ptb_start = ex.Trials(smp_idx).TrialStart_remappedGetSecs;
        end
        if isempty(ptb_start) || isnan(ptb_start)
            ptb_start = ex.Trials(smp_idx).TrialStart;
        end
        add_clock_align = ephys_start - ptb_start;
        [ex.Trials(smp_idx).Start] = [ex.Trials(smp_idx).Start] + add_clock_align;
        ex.Trials(smp_idx).TrialStart = ex.Trials(smp_idx).TrialStart + add_clock_align;
        ex.Trials(smp_idx).times.rewardGiven = ex.Trials(smp_idx).times.rewardGiven + add_clock_align;
        ex.Trials(smp_idx).times.reward = ex.Trials(smp_idx).times.reward + add_clock_align;
        S(iTrial).sample_onset = [S(iTrial).sample_onset ex.Trials(smp_idx).Start(1)];
        % eye data
        for nd = 1:length(ex.eyeCal.Delta)
            % get offset positions from online recalibration
            if smp_idx>=ex.eyeCal.Delta(nd).TrialNo
                RX0 = ex.eyeCal.Delta(nd).RX0;
                RY0 = ex.eyeCal.Delta(nd).RY0;
                LX0 = ex.eyeCal.Delta(nd).LX0;
                LY0 = ex.eyeCal.Delta(nd).LY0;
            end
        end
        % identify times where there is no eye-tracker data- either due
        % to blinks/closed eyes/lack of corneal reflection as the eye
        % moved outside the recording window. these are saved as very
        % high voltage values (±5)
        artefact_times = find((abs(tr(smp_idx).Eye.v(1,:)) > 4.9) | (abs(tr(smp_idx).Eye.v(2,:)) > 4.9) | (abs(tr(smp_idx).Eye.v(4,:)) > 4.9) | (abs(tr(smp_idx).Eye.v(5,:)) > 4.9));
        tr(smp_idx).Eye.v(1,artefact_times) = nan;
        tr(smp_idx).Eye.v(2,artefact_times) = nan;
        tr(smp_idx).Eye.v(3,artefact_times) = nan;
        tr(smp_idx).Eye.v(4,artefact_times) = nan;
        tr(smp_idx).Eye.v(5,artefact_times) = nan;
        tr(smp_idx).Eye.v(6,artefact_times) = nan;
        % convert voltages to DVA
        tr(smp_idx).Eye.v(1,:) = (tr(smp_idx).Eye.v(1,:)-ex.eyeCal.RX0 - RX0)*ex.eyeCal.RXGain*dpp;
        tr(smp_idx).Eye.v(2,:) = (tr(smp_idx).Eye.v(2,:)-ex.eyeCal.RY0 - RY0)*ex.eyeCal.RYGain*dpp;
        tr(smp_idx).Eye.v(4,:) = (tr(smp_idx).Eye.v(4,:)-ex.eyeCal.LX0 - LX0)*ex.eyeCal.LXGain*dpp;
        tr(smp_idx).Eye.v(5,:) = (tr(smp_idx).Eye.v(5,:)-ex.eyeCal.LY0 - LY0)*ex.eyeCal.LYGain*dpp;
        eye_x = (tr(smp_idx).Eye.v(1,:) + tr(smp_idx).Eye.v(4,:))*0.5; % mean x-position for both eyes
        eye_y = (tr(smp_idx).Eye.v(2,:) + tr(smp_idx).Eye.v(5,:))*0.5; % mean y-position for both eyes
        eye_xy_diff = [vecnorm(diff([eye_x; eye_y], 1, 2)) NaN]; % eye speed
        eye_r = sqrt(eye_x.^2 + eye_y.^2); % mean radial eye position for both eyes
        pupil = (tr(smp_idx).Eye.v(3,:) + tr(smp_idx).Eye.v(6,:))*0.5; % pupil size
        pupil_der = [diff(pupil) NaN]; % pupil derivative- marker for arousal. See Reimer et al. 2014
        eye_clock_align = ephys_start - ex.Trials(smp_idx).TrialStartDatapixx;
        eyetime = tr(smp_idx).Eye.t  + eye_clock_align; % time
        eyevals = [eyetime; eye_x; eye_y; eye_r; eye_xy_diff; pupil; pupil_der];
        S(iTrial).eyevals = [S(iTrial).eyevals eyevals];
    end
    % save the stimulus sequence values for each trial
    S(iTrial).cont_seq = [];
    for smp_idx = trl_smp_idx
        if strcmp(exname, 'M1_0057') % this session has an exception in how the values are stored
            if round(10000*ex.Trials(smp_idx).co)/10000 <= 0.0625
                ex.Trials(smp_idx).co = 0.0625;
            elseif round(10000*ex.Trials(smp_idx).co)/10000 > 0.0625 && round(10000*ex.Trials(smp_idx).co)/10000 <= 0.25
                ex.Trials(smp_idx).co = 0.25;
            elseif round(10000*ex.Trials(smp_idx).co)/10000 > 0.25 && round(10000*ex.Trials(smp_idx).co)/10000 <= 1
                ex.Trials(smp_idx).co = 1;
            elseif round(10000*ex.Trials(smp_idx).co)/10000 > 1
                ex.Trials(smp_idx).co = 1001;
            end
        end
        S(iTrial).cont_seq = [S(iTrial).cont_seq ex.Trials(smp_idx).co*ones(1, length(ex.Trials(smp_idx).Start))];
    end

    % get lags from fixation onset to stimulus onset, and inter-sample-interval for book-keeping
    S(iTrial).fix_lag = [(ex.Trials(trl_smp_idx(2)).Start(1)-ex.Trials(trl_smp_idx(2)).TrialStart) (ex.Trials(trl_smp_idx(3)).Start(1)-ex.Trials(trl_smp_idx(3)).TrialStart) (ex.Trials(trl_smp_idx(4)).Start(1)-ex.Trials(trl_smp_idx(4)).TrialStart)];
    S(iTrial).isi_lag = [(ex.Trials(trl_smp_idx(2)).TrialStart-ex.Trials(trl_smp_idx(1)).Start(end)) (ex.Trials(trl_smp_idx(3)).TrialStart-ex.Trials(trl_smp_idx(2)).Start(end)) (ex.Trials(trl_smp_idx(4)).TrialStart-ex.Trials(trl_smp_idx(3)).Start(end))];

    % save stimulus values and times
    S(iTrial).cont_time = [ex.Trials(trl_smp_idx).Start];
    S(iTrial).cont = [ex.Trials(trl_smp_idx).co];
    S(iTrial).lum = [ex.Trials(trl_smp_idx).st];

    S(iTrial).smp1 = [ex.Trials(trl_smp_idx(1)).Start];
    S(iTrial).smp2 = [ex.Trials(trl_smp_idx(2)).Start];
    S(iTrial).smp3 = [ex.Trials(trl_smp_idx(3)).Start];
    S(iTrial).smp4 = [ex.Trials(trl_smp_idx(4)).Start];

    S(iTrial).reward_size = ex.Trials(trl_num_idx).RewardSize;
    S(iTrial).reward_offset = ex.Trials(trl_num_idx).times.rewardGiven;
    S(iTrial).reward_onset = ex.Trials(trl_num_idx).times.reward;
    % get number of frames in each sample for book-keeping
    S(iTrial).nframes = [length([ex.Trials(trl_smp_idx(1)).Start]) length([ex.Trials(trl_smp_idx(2)).Start]) length([ex.Trials(trl_smp_idx(3)).Start]) length([ex.Trials(trl_smp_idx(4)).Start])];
    % --- timing
    S(iTrial).stim_start = S(iTrial).cont_time(1);

    % we need this variable in order to remove overlapped timepoints across
    % successive trials
    if trl_num_idx ~= size(ex.Trials, 2)
        nxt_trl_smp1 = trl_num_idx + 1;
        ephys_start = ex.Trials(nxt_trl_smp1).nevTrialStart;
        if strcmp(exname, 'M1_0037')
            ptb_start = ex.Trials(nxt_trl_smp1).TrialStart;
        else
            ptb_start = ex.Trials(nxt_trl_smp1).TrialStart_remappedGetSecs;
        end
        if isempty(ptb_start) || isnan(ptb_start)
            ptb_start = ex.Trials(nxt_trl_smp1).TrialStart;
        end
        add_clock_align = ephys_start - ptb_start;
        if ~isempty(ex.Trials(nxt_trl_smp1).Start)
            S(iTrial).nxt_trl_smp1 = [ex.Trials(nxt_trl_smp1).Start] + add_clock_align;
        else
            S(iTrial).nxt_trl_smp1 = ephys_start;
        end
    elseif trl_num_idx == size(ex.Trials, 2)
        S(iTrial).nxt_trl_smp1 = S(iTrial).smp4(end) + 3;
    end
end

fprintf('Binning and processing for regression analysis\n')

dt = ip.Results.dt;
num_pre = ceil(ip.Results.time_pre_stim/dt);
num_post = ceil(ip.Results.time_per_trial/dt);
% define the bins for each trial
frame_bins = -(num_pre-1)*dt:dt:(num_post-1)*dt;

% initialise the data structure
D = struct('cont_levels', unique(cell2mat(arrayfun(@(x) x.cont_seq, S', 'uni', 0))), ...
    'cont_seq', [], ...
    'smp1', [], ...
    'smp2', [], ...
    'smp3', [], ...
    'smp4', [], ...
    'isi_idx', [], ...
    'stim_timepts', [], ...
    'nonstim_timepts', [], ...
    'stim_start', [], ...
    'time_stamps', [], ...
    'trial_times', [], ...
    'reward_start', [], ...
    'robs', [], ...
    'blocks', [], ...
    'eye_pos', [], ...
    'eye_speed', [], ...
    'pupil', [], ...
    'pupil_derivative', []);

if exist('faceSVD', 'var')
    includeSVD = true;
    D.faceSVD = [];
    D.bodySVD = [];
else
    includeSVD = false;
end

% get unit probe, electrode, cluster id
pat = 'ch(?<probe>\w+)_(?<channel>\w+)_sort';
D.units = cellfun(@(x) regexp(x, pat, 'names'), fNm, 'uni', 1);
D.frameRate = 1/dt;
D.dt = dt;
binCtr = 0; % for keeping tracks of valid blocks

for iTrial = 1:nTrials

    fprintf('Trial %d/%d\n', iTrial, nTrials)

    % do some sanity checks and skip incomplete trials and exceptions
    % e.g., trials with too long fixation lag or inter-sample-interval

    if isempty(S(iTrial).cont_time)
        continue
    end
    if abs(S(iTrial).reward_size) == 0
        continue
    end
    if any(S(iTrial).fix_lag > 0.05)
        continue
    end
    if any(S(iTrial).isi_lag > 0.3)
        continue
    end
    if any(S(iTrial).nframes < 53)
        continue
    end

    % create bin times for this trial: we consider a non-overlapping bin
    % width specified by dt
    ft = [S(iTrial).smp1 S(iTrial).smp2 S(iTrial).smp3 S(iTrial).smp4];
    bin_times = ft(1) + frame_bins;
    assert(ft(end) < bin_times(end), 'Frames last longer than sample length! Use a longer sample length for sample')
    bin_edges = [bin_times bin_times(end)+dt];

    % get indices for this trial from the stored video SV timecourses
    meix = vid.t >= bin_times(1) & vid.t <= bin_times(end)+.1;
    t = vid.t(meix);
    % skip the trial if there is not video data available (just for
    % completeness)
    if isempty(t) || (t(1) - bin_times(1) > 0.3)
        continue
    end

    % identify the indices for inter-stimulus intervals
    isi_idx = [find(bin_times > S(iTrial).smp1(end) + 0.05 & bin_times < S(iTrial).smp2(1) - dt) ...
        find(bin_times > S(iTrial).smp2(end) + 0.05 & bin_times < S(iTrial).smp3(1) - dt) ...
        find(bin_times > S(iTrial).smp3(end) + 0.05 & bin_times < S(iTrial).smp4(1) - dt) ...
        find(bin_times > S(iTrial).smp4(end) + 1.05)];

    % identify the indices for stimulus presentation times: these
    % constitute our controlled retinal input epochs
    stim_timepts = [find(bin_times >= S(iTrial).smp1(1) - dt & bin_times <= S(iTrial).smp1(end) + 0.05) ...
        find(bin_times >= S(iTrial).smp2(1) - dt & bin_times <= S(iTrial).smp2(end) + 0.05) ...
        find(bin_times >= S(iTrial).smp3(1) - dt & bin_times <= S(iTrial).smp3(end) + 0.05) ...
        find(bin_times >= S(iTrial).smp4(1) - dt & bin_times <= S(iTrial).smp4(end) + 0.05)];

    % identify the indices for times in a trial that are not stimulus
    % epochs: these constitute the uncontrolled retinal input epochs
    nonstim_timepts = [find(bin_times < S(iTrial).smp1(1) - dt) ...
        find(bin_times > S(iTrial).smp4(end) + 0.05 & bin_times <= min(S(iTrial).smp4(end) + 1.05, S(iTrial).nxt_trl_smp1(1)))];

    % get the spike counts in each bin
    nbins = numel(bin_times);
    iix = st >= bin_times(1) & st <= bin_times(end);
    robs = zeros(nbins, NC, 'uint8');
    for cc = 1:NC
        robs(:,cc) = uint8(histcounts(st(iix & clu==cc), bin_edges));
    end

    % convert stim sequence to "One Hot" version
    cont_seq = double(S(iTrial).cont_seq(:)==D.cont_levels);
    [~, ~, cont_idx] = histcounts(ft, bin_edges);

    % re-bin at specified time-resolution
    cont_seq_hot = zeros(numel(bin_times), numel(D.cont_levels));
    cont_seq_hot(cont_idx, :) = cont_seq;

    % get the indices for onset of individual samples
    smp1 = histcounts(S(iTrial).smp1, bin_edges)';
    stimind = find(smp1);
    smp1(stimind) = 1;

    smp2 = histcounts(S(iTrial).smp2, bin_edges)';
    stimind = find(smp2);
    smp2(stimind) = 1;

    smp3 = histcounts(S(iTrial).smp3, bin_edges)';
    stimind = find(smp3);
    smp3(stimind) = 1;

    smp4 = histcounts(S(iTrial).smp4, bin_edges)';
    stimind = find(smp4);
    smp4(stimind) = 1;

    % get face and body SV timecourses
    meix = vid.t >= bin_times(1) & vid.t <= bin_times(end)+.1;
    t = vid.t(meix);
    if includeSVD
        faceSV = interp1(t, faceSVD.svdME(meix,:), bin_times, 'linear','extrap');
        D.faceSVD = [D.faceSVD; faceSV];
        bodySV = interp1(t, bodySVD.svdME(meix,:), bin_times, 'linear','extrap');
        D.bodySVD = [D.bodySVD; bodySV];
    end

    % get indices for onset of trial
    stim = histcounts(S(iTrial).stim_start, bin_edges)';
    stimind = find(stim);
    stim(stimind) = 1;

    % get indices for reward onset
    rwd = histcounts(S(iTrial).reward_onset, bin_edges)';
    rwdind = find(rwd);
    rwd(rwdind) = 1;

    % get eye data for this trial
    et = S(iTrial).eyevals(1, :); % time
    ex = S(iTrial).eyevals(2, :); % x position
    ey = S(iTrial).eyevals(3, :); % y position
    er = S(iTrial).eyevals(4, :); % radial position
    es = S(iTrial).eyevals(5, :); % speed
    ep = S(iTrial).eyevals(6, :); % pupil
    ep_der = S(iTrial).eyevals(7, :); % pupil derivative
    % get indices for this trial
    eix = et >= bin_times(1) & et <= bin_times(end)+.1;
    et2 = et(eix);
    eye_pos(1, :) = interp1(et2, ex(eix), bin_times, 'linear');
    eye_pos(2, :) = interp1(et2, ey(eix), bin_times, 'linear');
    eye_pos(3, :) = interp1(et2, er(eix), bin_times, 'linear');
    pupil = interp1(et2, ep(eix), bin_times, 'linear');
    eye_speed = interp1(et2, es(eix), bin_times, 'linear');
    pupil_derivative = interp1(et2, ep_der(eix), bin_times, 'linear');

    % update data structure
    D.time_stamps = [D.time_stamps; bin_times'];
    D.trial_times = [D.trial_times; bin_times(1) bin_times(end)];
    D.cont_seq = [D.cont_seq; uint8(cont_seq_hot)];
    D.robs = [D.robs; robs];
    D.smp1 = [D.smp1; uint8(smp1)];
    D.smp2 = [D.smp2; uint8(smp2)];
    D.smp3 = [D.smp3; uint8(smp3)];
    D.smp4 = [D.smp4; uint8(smp4)];
    D.stim_start = [D.stim_start; uint8(stim)];
    D.reward_start = [D.reward_start; uint8(rwd)];
    D.blocks = [D.blocks; [binCtr+1 binCtr+1+nbins]];
    D.isi_idx = [D.isi_idx; binCtr+isi_idx'];
    D.stim_timepts = [D.stim_timepts; binCtr+stim_timepts'];
    D.nonstim_timepts = [D.nonstim_timepts; binCtr+nonstim_timepts'];
    D.eye_pos = [D.eye_pos; eye_pos'];
    D.pupil = [D.pupil; pupil'];
    D.eye_speed = [D.eye_speed; eye_speed'];
    D.pupil_derivative = [D.pupil_derivative; pupil_derivative'];
    binCtr = binCtr + nbins;
end

% try to reduce memory footprint (runs on laptops)
D.cont_seq = uint8(D.cont_seq);
D.smp1 = uint8(D.smp1);
D.smp2 = uint8(D.smp2);
D.smp3 = uint8(D.smp3);
D.smp4 = uint8(D.smp4);
D.time_stamps = double(D.time_stamps);
D.trial_times = double(D.trial_times);
D.faceSVD = int16(D.faceSVD);
D.bodySVD = int16(D.bodySVD);
D.stim_start = uint8(D.stim_start);
D.reward_start = uint8(D.reward_start);
D.trial_time = frame_bins;
D.eye_pos_label = {'x', 'y', 'r'};
% also add receptive field data from Gabor fits
D.rf.x0 = RFposGab(:, 1);
D.rf.y0 = RFposGab(:, 2);
D.rf.r = sqrt(D.rf.x0.^2 + D.rf.y0.^2);
D.rf.width_x0 = RFposGab(:, 3);
D.rf.width_y0 = RFposGab(:, 4);
% add display related info
D.dpp = dpp;
D.screen.width = screen_width*dpp;
D.screen.height = screen_height*dpp;

% clean data
% this step is important because if we define a long trial length, like our
% default value of 5s, there is a danger that the post-stimulus interval of
% the current pseudo-trial will overlap with the stimulus onset or
% pre-stimulus interval of the next trial. To avoid fitting the same
% timepoints twice, we will identify the indices of 'repeated' timepoints
% and save them so we can remove them during model fitting
trial_times = D.trial_times;
time_stamps = D.time_stamps;
timepts2use = [];
for i = 1:size(D.blocks, 1)-1
    trl_end_idx = D.blocks(i, 2)-1;
    trl_start_idx = D.blocks(i, 1);
    nxt_trl_start_time = trial_times(i+1, 1);
    trl_end_time = trial_times(i, 2);
    if trl_end_time < nxt_trl_start_time
        new_trl_end_idx = trl_end_idx;
    else
        trl_time = time_stamps(trl_start_idx:trl_end_idx);
        new_trl_end_time = interp1(trl_time, trl_time, nxt_trl_start_time, 'previous');
        new_trl_end_idx = find(time_stamps == new_trl_end_time);
    end
    timepts2use = [timepts2use trl_start_idx:new_trl_end_idx];
end

% also add the last trial
trl_start_idx = D.blocks(i+1, 1);
trl_end_idx = D.blocks(i+1, 2)-1;
new_trl_end_idx = trl_end_idx;
timepts2use = [timepts2use trl_start_idx:new_trl_end_idx];
% remove the indices of variable inter-sample intervals
timepts2use = setdiff(timepts2use, D.isi_idx);
timepts2remove = setdiff(1:length(time_stamps), timepts2use);

% get the start and end indices of individual trials
D.blocks_clean = D.blocks;
D.blocks_clean(:, 2) = interp1(timepts2use, timepts2use, D.blocks(:, 2), 'previous', 'extrap');
D.blocks_clean(:, 1) = interp1(timepts2use, timepts2use, D.blocks(:, 1), 'next');
[~, D.blocks_clean(:, 1)] = intersect(timepts2use, D.blocks_clean(:, 1));
[~, D.blocks_clean(:, 2)] = intersect(timepts2use, D.blocks_clean(:, 2));
time_stamps = time_stamps(timepts2use);
% one last sanity check
assert(all(diff(time_stamps) > 0), 'time stamps was not fixed properly: there is a problem')
% get the start and end times of individual trials
D.trial_times_clean = D.trial_times;
D.trial_times_clean(:, 2) = time_stamps(D.blocks_clean(:, 2));
D.trial_times_clean(:, 1) = time_stamps(D.blocks_clean(:, 1));
% now save all the updated indices information
D.timepts2use = timepts2use;
D.timepts2remove = timepts2remove;
D.time_stamps = time_stamps;
D.nonstim_timepts = setdiff(D.nonstim_timepts, timepts2remove);
D.stim_timepts = setdiff(D.stim_timepts, timepts2remove);

% now save area information
V1_idx = find(unitArea == 1);
V2_idx = find(unitArea == 2);
V3_idx = find(unitArea == 3);
unassigned_idx = find(unitArea == 0);
if ~isempty(V1_idx)
    D.area.V1_idx = V1_idx;
end
if ~isempty(V2_idx)
    D.area.V2_idx = V2_idx;
end
if ~isempty(V3_idx)
    D.area.V3_idx = V3_idx;
end
if ~isempty(unassigned_idx)
    D.area.unassigned = unassigned_idx;
end
clear V1_idx V2_idx V3_idx unassigned_idx

% linearly translate receptive field positions according to the fixation
% and get a criteria for timepoints in each unit for which the translated
% rf positions are within the screen
% rf_criteria = 1: input to the RF can be inferred;
%                  RF within the boundaries of the display screen
% rf_criteria = 0: input to the RF cannot be inferred;
%                  eyes closed or RF outside the boundaries of the display
idx2use = find(~isnan(D.eye_pos(:, 3)));
D.rf_criteria = zeros(size(D.eye_pos, 1), size(D.robs, 2));
for i = 1:size(D.robs, 2)
    rel_rf_pos_x = D.rf.x0(i)+D.eye_pos(idx2use, 1);
    rel_rf_pos_y = D.rf.y0(i)+D.eye_pos(idx2use, 2);
    D.rf_criteria(idx2use, i) = (abs(rel_rf_pos_x) + 2 * D.rf.width_x0(i) < D.screen.width) & (abs(rel_rf_pos_y) + 2 * D.rf.width_y0(i) < D.screen.height);
end

% if we are including eye regressors, then we need to remove the
% timepoints where eye-data is not available before fitting
D_eye = D;
timepts2remove = union(find(isnan(D_eye.pupil) | isnan(D_eye.pupil_derivative)), D_eye.timepts2remove);
timepts2use = setdiff(D_eye.timepts2use, timepts2remove);
D_eye.blocks_clean = D_eye.blocks;
D_eye.blocks_clean(:, 2) = interp1(timepts2use, timepts2use, D_eye.blocks(:, 2), 'previous', 'extrap');
D_eye.blocks_clean(:, 1) = interp1(timepts2use, timepts2use, D_eye.blocks(:, 1), 'next', 'extrap');
[~, D_eye.blocks_clean(:, 1)] = intersect(timepts2use, D_eye.blocks_clean(:, 1));
[~, D_eye.blocks_clean(:, 2)] = intersect(timepts2use, D_eye.blocks_clean(:, 2));
D_eye.timepts2use = timepts2use;
D_eye.timepts2remove = timepts2remove;
D_eye.nonstim_timepts = setdiff(D_eye.nonstim_timepts, timepts2remove);
D_eye.stim_timepts = setdiff(D_eye.stim_timepts, timepts2remove);
end