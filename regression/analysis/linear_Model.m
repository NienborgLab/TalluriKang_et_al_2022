function linear_Model(exname, varargin)
% Bharath Talluri & Jacob Yates, Dec 2021
% adapted from linear encoding model used in the study "Single-trial neural
% dynamics are dominated by richly varied movements" by Musall et al., 2019

% Inputs:
%   datadir <char>: path to experiment director
%   exname  <char>: experiment name (e.g., ma_0823)

% Optional Inputs (as argument pairs)
%   'num_sv' <int>: number of video SVs to include in the analysis
%   'time_pre_stim'  <double>: time (sec) to include before stimulus onset
%   'video_time_embedding' <bool>: time (sec), duration of each trial
%   'include_eye' <bool>: whether the analysis requires eye regressors

ip = inputParser();
ip.addParameter('num_sv', 30)
ip.addParameter('time_pre_stim', .3)
ip.addParameter('video_time_embedding', 1)
ip.addParameter('include_eye', 0)
ip.parse(varargin{:});

include_eye = ip.Results.include_eye;
num_sv = ip.Results.num_sv;
video_time_embedding = ip.Results.video_time_embedding;
% load preprocessed dataset
if ~include_eye
    load(sprintf('../data/preprocessed/%s/preprocessed_data.mat', exname), 'D');
else
    load(sprintf('../data/preprocessed/%s/preprocessed_data.mat', exname), 'D_eye');
    D = D_eye;clear D_eye;
end
% make a directory to save the analysis files
results_dir = sprintf('../data/analysis/%s', exname);
if ~exist(results_dir, 'dir')
    mkdir(results_dir)
end

dt = D.dt;
%% asign some basic options for the model

% First check that the trials are properly accounted for
opts.frameRate = D.frameRate;
blocks = D.blocks;
frames = blocks(:,2) - blocks(:,1);
assert(all(frames==unique(frames)), 'makeDesignMatrix: there are more than one trial length. this is a problem')
opts.framesPerTrial = unique(frames);
opts.dt = dt;
opts.twin = 0.2;

% Set up parameters of the regression
% There are three event types when building the design matrix.
% Event type 1 will span the rest of the current trial.
% Type 2 will span frames according to sPostTime.
% Type 3 will span frames before the event
% according to mPreTime and frames after the event according to mPostTime.
opts.sPostTime = ceil(.25 * opts.frameRate);   % follow stim events for sPostStim in frames (used for eventType 2)
opts.mPreTime = ceil(0.5 * opts.frameRate);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
opts.mPostTime = ceil(0.5 * opts.frameRate);   % follow motor events for mPostStim in frames (used for eventType 3)
opts.folds = 10; %nr of folds for cross-validation

%% Task Events

% This is the important bit, you have to design what events feed into the
% design matrix. This will depend on what the stimulus and task are and how
% you want to represent their interactions. You'll notice here that I
% outline two possible ways to represent the effect of attention.
if strcmp(exname(1:2), 'M2')
    [taskEvents, taskEventType, taskLabels] = M2.getTaskVariabels(D, 'attention_mode', 'full', 'include_contrast', false);
    % remove the zero disparity regressors to avoid linear dependence
    % between regressors
    zero_labels = find(contains(taskLabels, {'left_0.00', 'right_0.00'}));
    use_labels_idx = setdiff(1:length(taskLabels), zero_labels);
    taskEvents = taskEvents(:, use_labels_idx);
    taskEventType = taskEventType(:, use_labels_idx);
    taskLabels = taskLabels(:, use_labels_idx);
    [nonTaskEvents, nonTaskEventType, nonTaskLabels] = M2.getNonTaskVariabels(D);
elseif strcmp(exname(1:2), 'M1')
    [taskEvents, taskEventType, taskLabels] = M1.getTaskVariabels(D, 'position_mode', 'full');
    [nonTaskEvents, nonTaskEventType, nonTaskLabels] = M1.getNonTaskVariabels(D);
end

disp("Building the design matrices")

% make design matrix
[taskR, taskIdx] = makeDesignMatrix(taskEvents, taskEventType, opts); %make design matrix for stim variables
[nonTaskR, nonTaskIdx] = makeDesignMatrix(nonTaskEvents, nonTaskEventType, opts); %make design matrix for non-stim variables

% remove repeated/invalid timepoints from the analysis: these timepoints
% were defined in the preprocessing
taskR(D.timepts2remove, :) = [];
nonTaskR(D.timepts2remove, :) = [];

% we define the number of top SVD components to take from the face and body
% SVDs - these can be different, but for now I am using equal number for
% both face and body
num_face_svd = num_sv;
num_body_svd = num_sv;
% the video features
face_svd = D.faceSVD(:, 1:num_face_svd);
body_svd = D.bodySVD(:, 1:num_body_svd);

% introduce time-embedding by convolving video regressors with basis
% functions
if video_time_embedding
    movement_regressors = lagged_video_regressors(face_svd, opts);
    pre_face_svdR = squeeze(movement_regressors(1, :, :));
    face_svdR = squeeze(movement_regressors(2, :, :));
    post_face_svdR = squeeze(movement_regressors(3, :, :));

    movement_regressors = lagged_video_regressors(body_svd, opts);
    pre_body_svdR = squeeze(movement_regressors(1, :, :));
    body_svdR = squeeze(movement_regressors(2, :, :));
    post_body_svdR = squeeze(movement_regressors(3, :, :));

    % remove repeated/invalid timepoints from the video SVs
    pre_face_svdR(D.timepts2remove, :) = [];
    face_svdR(D.timepts2remove, :) = [];
    post_face_svdR(D.timepts2remove, :) = [];
    pre_body_svdR(D.timepts2remove, :) = [];
    body_svdR(D.timepts2remove, :) = [];
    post_body_svdR(D.timepts2remove, :) = [];

    % zscore video SVs
    pre_face_svdR = zscore(double(pre_face_svdR));
    face_svdR = zscore(double(face_svdR));
    post_face_svdR = zscore(double(post_face_svdR));
    pre_body_svdR = zscore(double(pre_body_svdR));
    body_svdR = zscore(double(body_svdR));
    post_body_svdR = zscore(double(post_body_svdR));
else
    face_svdR = face_svd;
    body_svdR = body_svd;
    face_svdR(D.timepts2remove, :) = [];
    body_svdR(D.timepts2remove, :) = [];
    face_svdR = zscore(double(face_svdR));
    body_svdR = zscore(double(body_svdR));
end

if include_eye
    % include eye speed regressor as a proxy for eye movements
    eye_speed = D.eye_speed;
    eye_x = D.eye_pos(:, 1);
    eye_y = D.eye_pos(:, 2);
    eye_pupil = D.pupil;
    eye_pupilDerivative = D.pupil_derivative;
    % use a shorter window for basis functions for eye regressors- the eye
    % traces are not as smooth as video traces
    opts.eye_twin = 0.1;
    eye_speed_regressors = lagged_eye_regressors(eye_speed, opts, D.timepts2use(D.blocks_clean));
    pre_eye_sR = eye_speed_regressors(:, 1);
    eye_sR = eye_speed_regressors(:, 2);
    post_eye_sR = eye_speed_regressors(:, 3);

    pre_eye_sR(D.timepts2remove, :) = [];
    eye_sR(D.timepts2remove, :) = [];
    post_eye_sR(D.timepts2remove, :) = [];

    eye_x_regressors = lagged_eye_regressors(eye_x, opts, D.timepts2use(D.blocks_clean));
    pre_eye_xR = eye_x_regressors(:, 1);
    eye_xR = eye_x_regressors(:, 2);
    post_eye_xR = eye_x_regressors(:, 3);

    pre_eye_xR(D.timepts2remove, :) = [];
    eye_xR(D.timepts2remove, :) = [];
    post_eye_xR(D.timepts2remove, :) = [];

    eye_y_regressors = lagged_eye_regressors(eye_y, opts, D.timepts2use(D.blocks_clean));
    pre_eye_yR = eye_y_regressors(:, 1);
    eye_yR = eye_y_regressors(:, 2);
    post_eye_yR = eye_y_regressors(:, 3);

    pre_eye_yR(D.timepts2remove, :) = [];
    eye_yR(D.timepts2remove, :) = [];
    post_eye_yR(D.timepts2remove, :) = [];

    eye_pupil_regressors = lagged_eye_regressors(eye_pupil, opts, D.timepts2use(D.blocks_clean));
    pre_eye_pR = eye_pupil_regressors(:, 1);
    eye_pR = eye_pupil_regressors(:, 2);
    post_eye_pR = eye_pupil_regressors(:, 3);

    pre_eye_pR(D.timepts2remove, :) = [];
    eye_pR(D.timepts2remove, :) = [];
    post_eye_pR(D.timepts2remove, :) = [];

    eye_pupilDer_regressors = lagged_eye_regressors(eye_pupilDerivative, opts, D.timepts2use(D.blocks_clean));
    pre_eye_pdR = eye_pupilDer_regressors(:, 1);
    eye_pdR = eye_pupilDer_regressors(:, 2);
    post_eye_pdR = eye_pupilDer_regressors(:, 3);

    pre_eye_pdR(D.timepts2remove, :) = [];
    eye_pdR(D.timepts2remove, :) = [];
    post_eye_pdR(D.timepts2remove, :) = [];

    % zscore
    pre_eye_sR = zscore(double(pre_eye_sR));
    eye_sR = zscore(double(eye_sR));
    post_eye_sR = zscore(double(post_eye_sR));
    pre_eye_xR = zscore(double(pre_eye_xR));
    eye_xR = zscore(double(eye_xR));
    post_eye_xR = zscore(double(post_eye_xR));
    pre_eye_yR = zscore(double(pre_eye_yR));
    eye_yR = zscore(double(eye_yR));
    post_eye_yR = zscore(double(post_eye_yR));
    pre_eye_pR = zscore(double(pre_eye_pR));
    eye_pR = zscore(double(eye_pR));
    post_eye_pR = zscore(double(post_eye_pR));
    pre_eye_pdR = zscore(double(pre_eye_pdR));
    eye_pdR = zscore(double(eye_pdR));
    post_eye_pdR = zscore(double(post_eye_pdR));
end

% get drift regressors
NT = size(D.robs,1);
% define the number of basis functions to include
if strcmp(exname(1:2), 'M1')
    nBasis = 10;
elseif strcmp(exname(1:2), 'M2')
    nBasis = 8;
end
step = ceil(NT/nBasis);
basis_centers = 0:step:NT;
sess_driftR = max(1 - abs( (1:NT)'-basis_centers)/step, 0);
driftR = sess_driftR;
driftR(D.timepts2remove, :) = [];

% get firing rates
Vc = double(D.robs'); % spike counts

% get indices for cross validation
% crossvalidation will be done at a single trial level
total_time = size(Vc,2); % total time points
dataIdx = true(opts.folds, total_time);
n = unique(frames); % bins per trial
trl_order = randperm(total_time / n); % randomly shuffle trials
for t = 1:(total_time/n) % loop over trials
    i = mod(t, opts.folds); % assign trial to fold
    if i == 0, i = opts.folds; end
    dataIdx(i,(trl_order(t)-1)*n + (1:n)) = false; % set this trial to test set
end

% remove repeated/invalid timepoints
Vc(:, D.timepts2remove) = [];
dataIdx(:, D.timepts2remove) = [];

if video_time_embedding
    if ~include_eye
        fullR = [taskR, nonTaskR, pre_face_svdR, face_svdR, post_face_svdR, pre_body_svdR, body_svdR, post_body_svdR]; %make new, single design matrix
        addLabels = [nonTaskLabels, {'pre_face_svd'}, {'face_svd'}, {'post_face_svd'}, {'pre_body_svd'}, {'body_svd'}, {'post_body_svd'}];
        regIdx = [taskIdx; nonTaskIdx + max(taskIdx)]; % add non-task indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(pre_face_svdR, 2), 1)]; % add pre_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(face_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(post_face_svdR, 2), 1)]; % add post_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(pre_body_svdR, 2), 1)]; % add pre_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(body_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(post_body_svdR, 2), 1)]; % add post_movement SVD indices
    else
        fullR = [taskR, nonTaskR, pre_face_svdR, face_svdR, post_face_svdR, pre_body_svdR, body_svdR, post_body_svdR, pre_eye_sR, eye_sR, post_eye_sR, pre_eye_xR, eye_xR, post_eye_xR, pre_eye_yR, eye_yR, post_eye_yR, pre_eye_pR, eye_pR, post_eye_pR, pre_eye_pdR, eye_pdR, post_eye_pdR]; %make new, single design matrix
        addLabels = [nonTaskLabels, {'pre_face_svd'}, {'face_svd'}, {'post_face_svd'}, {'pre_body_svd'}, {'body_svd'}, {'post_body_svd'}, {'pre_eyeSpeed'}, {'eyeSpeed'}, {'post_eyeSpeed'}, {'pre_eyePos_x'}, {'eyePos_x'}, {'post_eyePos_x'}, {'pre_eyePos_y'}, {'eyePos_y'}, {'post_eyePos_y'}, {'pre_pupil_size'}, {'pupil_size'}, {'post_pupil_size'}, {'pre_pupil_der'}, {'pupil_der'}, {'post_pupil_der'}];
        regIdx = [taskIdx; nonTaskIdx + max(taskIdx)]; % add non-task indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(pre_face_svdR, 2), 1)]; % add pre_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(face_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(post_face_svdR, 2), 1)]; % add post_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(pre_body_svdR, 2), 1)]; % add pre_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(body_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(post_body_svdR, 2), 1)]; % add post_movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_pupil derivative indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pupil derivative indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_pupil derivative indices
    end
else
    if ~include_eye
        fullR = [taskR, nonTaskR, face_svdR, body_svdR, pre_eye_sR, eye_sR, post_eye_sR, pre_eye_xR, eye_xR, post_eye_xR, pre_eye_yR, eye_yR, post_eye_yR, pre_eye_pR, eye_pR, post_eye_pR, pre_eye_pdR, eye_pdR, post_eye_pdR]; %make new, single design matrix
        addLabels = [nonTaskLabels, {'face_svd'}, {'body_svd'}, {'pre_eyeSpeed'}, {'eyeSpeed'}, {'post_eyeSpeed'}, {'pre_eyePos_x'}, {'eyePos_x'}, {'post_eyePos_x'}, {'pre_eyePos_y'}, {'eyePos_y'}, {'post_eyePos_y'}, {'pre_pupil_size'}, {'pupil_size'}, {'post_pupil_size'}, {'pre_pupil_der'}, {'pupil_der'}, {'post_pupil_der'}];
        regIdx = [taskIdx; nonTaskIdx + max(taskIdx)]; % add non-task indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(face_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(body_svdR, 2), 1)]; % add movement SVD indices
    else
        fullR = [taskR, nonTaskR, face_svdR, body_svdR]; %make new, single design matrix
        addLabels = [nonTaskLabels, {'face_svd'}, {'body_svd'}];
        regIdx = [taskIdx; nonTaskIdx + max(taskIdx)]; % add non-task indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(face_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, size(body_svdR, 2), 1)]; % add movement SVD indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye speed indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye x indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_eye y indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_pupil indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pre_pupil derivative indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add pupil derivative indices
        regIdx = [regIdx; repmat(max(regIdx)+1, 1, 1)]; % add post_pupil derivative indices
    end
end
regLabels = [taskLabels, addLabels];

% remove empty regressors
cIdx = find(sum(abs(fullR),1) == 0); %don't use empty regressors
fullR(:, cIdx) = [];
regIdx(cIdx) = [];

% now add drift
fullR = [fullR, driftR]; %make new, single design matrix
regLabels = [regLabels, {'drift'}];
regIdx = [regIdx; repmat(max(regIdx)+1, size(driftR, 2), 1)]; % add drift indices

% to get the task model, we shuffle the video svd regressors, to keep the
% dimensionality constant
shuffle_regressors = 'svd';
shuffle_idx = find(ismember(regIdx, find(contains(regLabels, shuffle_regressors))));
ModR = fullR;
for idx = shuffle_idx'
    ModR(:, idx) = datasample(ModR(:, idx), size(ModR, 1), 1);
end

% save variables
if ~include_eye
    save(sprintf('%s/analysisVars.mat', results_dir), 'Vc', 'opts', 'fullR', 'regIdx', 'regLabels', 'ModR', 'dataIdx',  '-v7.3');
else
    save(sprintf('%s/analysisVars_withEye.mat', results_dir), 'Vc', 'opts', 'fullR', 'regIdx', 'regLabels', 'ModR', 'dataIdx',  '-v7.3');
end
disp("Done")

%% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% fully orthogonal to all preceeding regressors in the matrix and 0 being
% fully redundant. Having fully redundant regressors in the matrix will
% break the model, so in this example those regressors are removed. In
% practice, you should understand where the redundancy is coming from and
% change your model design to avoid it in the first place!

[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
diagQRR = abs(diag(fullQRR));

figure(1); clf
set(gcf, 'DefaultAxesColorOrder', hsv(max(regIdx)));
for i = unique(regIdx)'
    inds = find(regIdx==i);
    h(i) = plot(inds, diagQRR(inds),'linewidth',2); hold on;
    ylim([-0.1 1.1]);
    title('Regressor orthogonality');
end
legend(h, regLabels, 'Location', 'BestOutside', 'FontSize', 5);
drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; box("off");
ylabel('Norm. vector angle');
xlabel('Design matrix columns');
assert(sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) >= size(fullR,2), 'Design matrix is rank-defficient. Check the design matrix for linear dependence between columns.\n');
set(gcf, 'PaperType', 'A4', 'PaperOrientation','landscape');
print(gcf, '-dpdf', '-fillpage', sprintf('%s/regressor_orthogonality.pdf', results_dir));

%% fit linear model and compute variance explained
if strcmp(exname(1:2), 'M1')
    fit_linear_Model_M1(exname, include_eye);
elseif strcmp(exname(1:2), 'M2')
    fit_linear_Model_M2(exname, include_eye);
end
