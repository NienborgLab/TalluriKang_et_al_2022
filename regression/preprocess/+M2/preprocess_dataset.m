function [D, D_eye] = preprocess_dataset(datadir, exname, varargin)
% first written by J. Yates, 2021
% modified by Bharath Talluri, Dec 2021
% Modified for sharing by Bharath Talluri, Dec 2022

% Preprocess dataset will take all of the files in an experiment folder,
% time-align them for each trial and then bin and make it time-continuous

% Inputs:
%   datadir <char>: path to experiment director
%   exname  <char>: experiment name (e.g., ma_0823)

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
%   faceSVD        <int16>: NT x 200, SVD components of face motion energy
%   bodySVD        <int16>: NT x 200, SVD components of body motion energy
%   units         <strcut>: struct array of unit information
%   frameRate     <double>: time resolution


ip = inputParser();
ip.addParameter('dt', 1/60)
ip.addParameter('time_pre_stim', .3)
ip.addParameter('time_per_trial', 3)
ip.parse(varargin{:});

% load data
fl = dir(fullfile(datadir, [exname '*']));
fpath = fullfile(fl.folder, fl.name);

flist = {dir(fullfile(fpath, '*.mat')).name}; % changed this to make sure this runs on windows
for f = 1:numel(flist)
    if ~isempty(flist{f})
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
            load(fullfile(fpath, flist{f})) %#ok<*LOAD>
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

nTrials = numel(ex.Trials);

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

% S is a struct array of length N Trials and accomplishes time-aligning
% everything within each trial

S = repmat(struct('hdx', [], ...
    'hdx2', [], ...
    'hdx_time', [], ...
    'hdx_seq', [], ...
    'hdx_seq2', [], ...
    'dc', [], ...
    'dc2', [], ...
    'targ_good', [], ...
    'xcue', [], ...
    'reward', [], ...
    'resp', [], ...
    'targ_time', [], ...
    'choice', [], ...
    'choice_time', [], ...
    'correct', [], ...
    'stim_onset_time', [], ...
    'reward_onset_time', [], ...
    'fp_onset', [], ...
    'fp_offset', [], ...
    'fix_start', [], ...
    'fix_complete', [], ...
    'break_fix', [], ...
    'frame_times', [], ...
    'eyevals', []), nTrials, 1);

fprintf('Get stimulus information \n');

for iTrial = 1:nTrials
    fprintf('Trial %d/%d\n', iTrial, nTrials)

    % times all frames flipped according to psychtoolbox
    frame_times = ex.Trials(iTrial).flip_info(:,2);

    % clock alignment: there are different clocks in the experimental
    % session: the clock from stimulus computer, the clock on the neural recording
    % computer, and the clock on the eye-tracker computer. We need to make
    % sure these clocks are aligned before we can align spikes, stimulus
    % events, and eye/pupil timeseries
    % align the clocks and save the sample onset times in the data
    % structure
    ephys_start = ex.Trials(iTrial).nevTrialStart;
    ptb_start = ex.Trials(iTrial).TrialStart_remappedGetSecs;
    tfun = @(x) x - ptb_start + ephys_start;

    % stimulus values for this trial
    S(iTrial).hdx_time = tfun(ex.Trials(iTrial).Start);
    S(iTrial).hdx = ex.Trials(iTrial).hdx;
    S(iTrial).hdx2 = ex.Trials(iTrial).hdx2;
    S(iTrial).hdx_seq = round(ex.Trials(iTrial).hdx_seq(1:ex.Trials(iTrial).framecnt)*1000)/1000;
    if isnan(ex.Trials(iTrial).hdx_seq2)
        S(iTrial).hdx_seq2 = nan(1, ex.Trials(iTrial).framecnt);
    else
        S(iTrial).hdx_seq2 = round(ex.Trials(iTrial).hdx_seq2(1:ex.Trials(iTrial).framecnt)*1000)/1000;
    end
    S(iTrial).dc = ex.Trials(iTrial).Dc;
    S(iTrial).dc2 = ex.Trials(iTrial).Dc2;
    S(iTrial).xcue = ex.Trials(iTrial).x0;
    S(iTrial).reward = ex.Trials(iTrial).Reward;
    S(iTrial).resp = ex.Trials(iTrial).RespDir;
    S(iTrial).targ_good = ex.Trials(iTrial).targ.goodT;
    S(iTrial).instruct = ex.Trials(iTrial).instructionTrial;
    if S(iTrial).reward == 1
        S(iTrial).choice = sign(S(iTrial).hdx);
    elseif S(iTrial).reward == -1
        S(iTrial).choice = -sign(S(iTrial).hdx);
    elseif S(iTrial).reward == 0
        S(iTrial).choice = NaN;
    end

    % timing information
    S(iTrial).frame_times = tfun(frame_times);
    % stim
    stimFrameOn = find( (frame_times - ex.Trials(iTrial).times.stimOn) > 0, 1, 'first');
    S(iTrial).stim_onset_time = tfun(frame_times(stimFrameOn));
    % reward
    S(iTrial).reward_onset_time = tfun(ex.Trials(iTrial).times.reward);
    % target
    targetFrameOn = find( (frame_times - ex.Trials(iTrial).times.targOn) > 0, 1, 'first');
    S(iTrial).target_onset = tfun(frame_times(targetFrameOn));
    % fixation
    fixOnset = find( (frame_times - ex.Trials(iTrial).times.fpOn) > 0, 1, 'first');
    S(iTrial).fp_onset = tfun(frame_times(fixOnset));
    S(iTrial).fix_start = tfun(ex.Trials(iTrial).times.startFixation);
    S(iTrial).break_fix = tfun(ex.Trials(iTrial).times.breakFixation);
    fixOffset = find( (frame_times - ex.Trials(iTrial).times.fpOff) > 0, 1, 'first');
    S(iTrial).fp_offset = tfun(frame_times(fixOffset));
    if isempty(S(iTrial).fp_offset) % use breakfix
        fixOffset = find( (frame_times - ex.Trials(iTrial).times.breakFixation) > 0, 1, 'first');
        S(iTrial).fp_offset = tfun(frame_times(fixOffset));
    end
    % choice
    S(iTrial).choice_time = tfun(ex.Trials(iTrial).times.choice);
    % target
    S(iTrial).targ_time = tfun(ex.Trials(iTrial).times.targOn);
    % stim onset
    if ~isempty(S(iTrial).hdx_time)
        S(iTrial).stim_start = S(iTrial).hdx_time(1);
    else
        S(iTrial).stim_start = [];
    end

    timingFields = {'stim_onset_time', 'reward_onset_time', 'target_onset', 'fp_onset', 'fp_offset', 'fix_start', 'fix_complete', 'break_fix', 'choice_time', 'targ_time'};
    for ifield = 1:numel(timingFields)
        if isempty(S(iTrial).(timingFields{ifield}))
            S(iTrial).(timingFields{ifield}) = nan;
        end
    end
    % eye data
    for nd = 1:length(ex.eyeCal.Delta)
        if iTrial>=ex.eyeCal.Delta(nd).TrialNo
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
    artefact_times = find((abs(tr(iTrial).Eye.v(1,:)) > 4.9) | (abs(tr(iTrial).Eye.v(2,:)) > 4.9) | (abs(tr(iTrial).Eye.v(4,:)) > 4.9) | (abs(tr(iTrial).Eye.v(5,:)) > 4.9));
    tr(iTrial).Eye.v(1,artefact_times) = nan;
    tr(iTrial).Eye.v(2,artefact_times) = nan;
    tr(iTrial).Eye.v(3,artefact_times) = nan;
    tr(iTrial).Eye.v(4,artefact_times) = nan;
    tr(iTrial).Eye.v(5,artefact_times) = nan;
    tr(iTrial).Eye.v(6,artefact_times) = nan;
    % convert voltages to DVA
    tr(iTrial).Eye.v(1,:) = (tr(iTrial).Eye.v(1,:)-ex.eyeCal.RX0 - RX0)*ex.eyeCal.RXGain*dpp;
    tr(iTrial).Eye.v(2,:) = (tr(iTrial).Eye.v(2,:)-ex.eyeCal.RY0 - RY0)*ex.eyeCal.RYGain*dpp;
    tr(iTrial).Eye.v(4,:) = (tr(iTrial).Eye.v(4,:)-ex.eyeCal.LX0 - LX0)*ex.eyeCal.LXGain*dpp;
    tr(iTrial).Eye.v(5,:) = (tr(iTrial).Eye.v(5,:)-ex.eyeCal.LY0 - LY0)*ex.eyeCal.LYGain*dpp;
    eye_x = (tr(iTrial).Eye.v(1,:) + tr(iTrial).Eye.v(4,:))*0.5; % mean x-position for both eyes
    eye_y = (tr(iTrial).Eye.v(2,:) + tr(iTrial).Eye.v(5,:))*0.5; % mean y-position for both eyes
    eye_xy_diff = [vecnorm(diff([eye_x; eye_y], 1, 2)) NaN]; % eye speed
    eye_r = sqrt(eye_x.^2 + eye_y.^2); % mean radial eye position for both eyes
    pupil = (tr(iTrial).Eye.v(3,:) + tr(iTrial).Eye.v(6,:))*0.5; % pupil size
    pupil_der = [diff(pupil) NaN]; % pupil derivative- marker for arousal. See Reimer et al. 2014
    eye_clock_align = ephys_start - ex.Trials(iTrial).TrialStartDatapixx;
    eyetime = tr(iTrial).Eye.t  + eye_clock_align;
    eyevals = [eyetime; eye_x; eye_y; eye_r; eye_xy_diff; pupil; pupil_der];
    S(iTrial).eyevals = [S(iTrial).eyevals eyevals];

    % we need this variable in order to remove overlapped timepoints across
    % successive trials
    if iTrial ~= size(ex.Trials, 2)
        nxt_iTrial = iTrial + 1;
        ephys_start = ex.Trials(nxt_iTrial).nevTrialStart;
        ptb_start = ex.Trials(nxt_iTrial).TrialStart_remappedGetSecs;
        add_clock_align = ephys_start - ptb_start;
        if ~isempty(ex.Trials(nxt_iTrial).Start)
            S(iTrial).nxt_iTrial = [ex.Trials(nxt_iTrial).Start] + add_clock_align;
        else
            S(iTrial).nxt_iTrial = ephys_start;
        end
    elseif iTrial == size(ex.Trials, 2)
        if ~isempty(S(iTrial).hdx_time)
            S(iTrial).nxt_iTrial = S(iTrial).hdx_time(end) + 3;
        else
            S(iTrial).nxt_iTrial = [];
        end
    end
end

fprintf('Binning and processing for regression analysis\n')

dt = ip.Results.dt;
num_pre = ceil(ip.Results.time_pre_stim/dt);
num_post = ceil(ip.Results.time_per_trial/dt);
% define the bins for each trial
frame_bins = -(num_pre-1)*dt:dt:(num_post-1)*dt;

% initialise the data structure
D = struct('hdx_levels', unique(cell2mat(arrayfun(@(x) x.hdx_seq, S', 'uni', 0))), ...
    'hdx_right', [], ...
    'hdx_left', [], ...
    'contrast_left', [], ...
    'contrast_right', [], ...
    'attend', [], ...
    'stim_start', [], ...
    'stim_timepts', [], ...
    'nonstim_timepts', [], ...
    'time_stamps', [], ...
    'trial_times', [], ...
    'stim_onset_time', [], ...
    'reward_onset_time', [], ...
    'reward_start', [], ...
    'robs', [], ...
    'blocks', [], ...
    'resp', [], ...
    'choice', [], ...
    'target_left', [], ...
    'target_right', [], ...
    'correct', [], ...
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
pat = '(?<probe>\w+).elec(?<electrode>\w+)\s+spk_(?<cluster>\w+)';
D.units = cellfun(@(x) regexp(x, pat, 'names'), fNm, 'uni', 1);
D.frameRate = 1/dt;
D.dt = dt;
binCtr = 0; % for keeping tracks of valid blocks

% check which hemisphere the units are in
D.units(1).rf_loc = {};

for i = 1:length(D.units)
    port = D.units(i).probe;
if ~strcmp(exname, 'M2_0824')
    probe_id = find(strcmp({ex.setup.gv.probe.port}, port));
    if ex.setup.gv.probe(probe_id).leftHemisphereRecorded
        D.units(i).rf_loc = 'right';
    elseif ex.setup.gv.probe(probe_id).rightHemisphereRecorded
        D.units(i).rf_loc = 'left';
    end
else
    if ex.setup.gv.probe.leftHemisphereRecorded
        D.units(i).rf_loc = 'right';
    elseif ex.setup.gv.probe.rightHemisphereRecorded
        D.units(i).rf_loc = 'left';
    end
end
end

for iTrial = 1:nTrials

    fprintf('Trial %d/%d\n', iTrial, nTrials)

    % do some sanity checks and skip incomplete trials and exceptions
    % e.g., instruction trials or trials that were not completed
    if isnan(S(iTrial).resp) || S(iTrial).instruct == 1
        continue
    end

    % create bin times for this trial: we consider a non-overlapping bin
    % width specified by dt
    ft = S(iTrial).hdx_time;
    bin_times = ft(1) + frame_bins;
    assert(ft(end) < bin_times(end), 'Frames last longer than trial length! Use a longer trial length')

    % get indices for this trial from the stored video SV timecourses
    meix = vid.t >= bin_times(1) & vid.t <= bin_times(end)+.1;
    t = vid.t(meix);
    % skip the trial if there is not video data available (just for
    % completeness)
    if isempty(t) || (t(1) - bin_times(1) > 0.3)
        continue
    end

    % identify the indices for stimulus presentation times: these
    % constitute our controlled retinal input epochs
    stim_timepts = find(bin_times >= S(iTrial).hdx_time(1) -dt & bin_times <= S(iTrial).hdx_time(1) + 2);

    % identify the indices for times in a trial that are not stimulus
    % epochs: these constitute the uncontrolled retinal input epochs
    nonstim_timepts = [find(bin_times < S(iTrial).hdx_time(1) - dt) ...
        find(bin_times > S(iTrial).hdx_time(1) + 2 & bin_times <= min(S(iTrial).hdx_time(1) + 3, S(iTrial).nxt_iTrial(1)))];

    % get the spike counts in each bin
    nbins = numel(bin_times);
    iix = st >= bin_times(1) & st <= bin_times(end);
    robs = zeros(nbins, NC, 'uint8');
    bin_edges = [bin_times bin_times(end)+dt];
    for cc = 1:NC
        robs(:,cc) = uint8(histcounts(st(iix & clu==cc), bin_edges));
    end

    % convert hdx sequence to "One Hot" version
    hdx_seq1 = double(S(iTrial).hdx_seq(:)==D.hdx_levels);
    hdx_seq2 = double(S(iTrial).hdx_seq2(:)==D.hdx_levels);
    [~, ~, hdx_idx] = histcounts(ft, bin_times);

    % re-bin at specified time-resolution
    hdx_seq_1hot = zeros(numel(bin_times), numel(D.hdx_levels));
    hdx_seq2_1hot = zeros(numel(bin_times), numel(D.hdx_levels));
    hdx_seq_1hot(hdx_idx, :) = hdx_seq1;
    hdx_seq2_1hot(hdx_idx, :) = hdx_seq2;

    % correctly account for the cue
    if S(iTrial).xcue < 0
        hdx_left = hdx_seq_1hot;
        hdx_right = hdx_seq2_1hot;
    else
        hdx_left = hdx_seq2_1hot;
        hdx_right = hdx_seq_1hot;
    end

    % keep track of when the stimulus was on
    contrast_left = sum(hdx_left, 2)>0;
    contrast_right = sum(hdx_right, 2)>0;

    % get face and body SV timecourses
    meix = vid.t >= bin_times(1) & vid.t <= bin_times(end)+.1;
    t = vid.t(meix);
    if includeSVD
        faceSV = interp1(t, faceSVD.svdME(meix,:), bin_times, 'linear','extrap');
        D.faceSVD = [D.faceSVD; faceSV];
        bodySV = interp1(t, bodySVD.svdME(meix,:), bin_times, 'linear','extrap');
        D.bodySVD = [D.bodySVD; bodySV];
    end

    % response
    cho = histcounts(S(iTrial).choice_time, bin_edges)';
    choind = find(cho);
    cho(choind) = sign(S(iTrial).resp-1.5);

    % binary choice
    bin_cho = histcounts(S(iTrial).choice_time, bin_edges)';
    bin_cho(choind) = S(iTrial).choice;

    % choice accuracy
    correct = cho;
    if ~isempty(S(iTrial).correct)
        correct(choind) = S(iTrial).correct;
    else
        correct = zeros(size(correct));
    end

    % targets
    targ = histcounts(S(iTrial).targ_time, bin_edges)';
    target_left = targ;
    target_right = targ;
    targind = find(targ);
    if S(iTrial).xcue < 0  % attend left
        target_left(targind) = 0;
        target_right(targind) = 1;
    else  % attend right
        target_left(targind) = 1;
        target_right(targind) = 0;
    end

    % stim start
    stim = histcounts(S(iTrial).stim_start, bin_edges)';
    stimind = find(stim);
    stim(stimind) = 1;

    % reward start
    rwd = histcounts(S(iTrial).reward_onset_time, bin_edges)';
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
    D.eye_pos = [D.eye_pos; eye_pos'];
    D.pupil = [D.pupil; pupil'];
    D.eye_speed = [D.eye_speed; eye_speed'];
    D.pupil_derivative = [D.pupil_derivative; pupil_derivative'];

    % update
    D.time_stamps = [D.time_stamps; bin_times'];
    D.trial_times = [D.trial_times; bin_times(1) bin_times(end)];
    D.hdx_right = [D.hdx_right; hdx_right];
    D.hdx_left = [D.hdx_left; hdx_left];
    D.contrast_left = [D.contrast_left; contrast_left];
    D.contrast_right = [D.contrast_right; contrast_right];
    D.robs = [D.robs; robs];
    D.resp = [D.resp; cho];
    D.choice = [D.choice; bin_cho];
    D.stim_start = [D.stim_start; uint8(stim)];
    D.reward_start = [D.reward_start; uint8(rwd)];
    D.correct = [D.correct; correct];
    D.target_left = [D.target_left; target_left];
    D.target_right = [D.target_right; target_right];
    D.attend = [D.attend; sign(S(iTrial).xcue)*ones(size(correct))];
    D.stim_timepts = [D.stim_timepts; binCtr+stim_timepts'];
    D.nonstim_timepts = [D.nonstim_timepts; binCtr+nonstim_timepts'];
    D.stim_onset_time = [D.stim_onset_time S(iTrial).stim_start];
    D.reward_onset_time = [D.reward_onset_time S(iTrial).stim_start + (rwdind - stimind)*dt];
    D.blocks = [D.blocks; [binCtr+1 binCtr+1+nbins]];
    binCtr = binCtr + nbins;
end

% try to reduce memory footprint (runs on laptops)
D.hdx_right = uint8(D.hdx_right);
D.hdx_left = uint8(D.hdx_left);
D.faceSVD = int16(D.faceSVD);
D.bodySVD = int16(D.bodySVD);
D.resp = int8(D.resp);
D.choice = int8(D.choice);
D.contrast_left = uint8(D.contrast_left);
D.stim_start = uint8(D.stim_start);
D.reward_start = uint8(D.reward_start);
D.contrast_right = uint8(D.contrast_right);
D.target_left = uint8(D.target_left);
D.target_right = uint8(D.target_right);
D.trial_time = frame_bins;
D.eye_pos_label = {'x', 'y', 'r'};
% also add receptive field data from Gabor fits
% receptive field data
D.rf.x0 = rfGabor.rfX;
D.rf.y0 = rfGabor.rfY;
D.rf.r = sqrt(D.rf.x0.^2 + D.rf.y0.^2);
D.rf.width_x0 = rfGabor.widthX;
D.rf.width_y0 = rfGabor.widthY;
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
timepts2remove = union(find(isnan(D_eye.pupil)), D_eye.timepts2remove);
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