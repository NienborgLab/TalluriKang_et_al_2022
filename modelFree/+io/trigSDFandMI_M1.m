function trigSDFandMI_M1(sesID,tsDir,trigDir,algTrigDir,rrmDir,analDir) 

% function trigSDFandMI_M1(sesID,tsDir,trigDir,algTrigDir,rrmDir,analDir) 
%
% input:    
% sesID         session ID
% tsDir         path to dir of preprocessed data for RRM
% trigDir       path to dir of triggerTimes for (no-)movement onset epochs
% algTrigDir    path to dir of indices to trigger times aligned to
%               time-stamps of preprocessed data
% rrmDir        path to dir of linear model results (for de-trending)
% analDir       path to dir to store results from current analysis 



% parameters
fS = 60;  % sampling frequency of data/model (60Hz)
tWin = [-1 1]; % timewindow for SDF (in sec; +/-1sec)
w = abs(tWin*fS); % width window around trigger of SDF
tWinsc = [0.15 0.45]; % timewindow for spike-counts to compute MotionIndex (MI)


fprintf('processing... %s \n',sesID)

% load trigger indices
load(fullfile(algTrigDir,sesID,'both_trigIdx'),...
    'trigIdxS','trigIdxSNM','trigIdxM','trigIdxNM');

% load trigger times
load(fullfile(trigDir,sesID,'both_TrigTimes'),'StM','StNMM');

% load data and model predictions 
load(fullfile(rrmDir,sesID,'analysisVars'),'Vc');
load(fullfile(rrmDir,sesID,'crossval_modelFits'),'Vfull','Vfull_drift');

% load time stamps
load(fullfile(tsDir,sesID,'preprocessed_data'),'D');    
ts = D.time_stamps;

V = Vc - Vfull_drift + mean(Vfull_drift,2); % drift corrected data
Vfl = Vfull - Vfull_drift + mean(Vfull_drift,2); % full model, drift corrected

nU = size(Vc,1); %# units
nS = length(trigIdxS); %#stimuli

% compute triggered SDFs ==================================================
sdfStM = cell(nU,nS); sdfStNMM = cell(nU,nS);  % neural data
sdfMM = cell(nU,nS); sdfNMM = cell(nU,nS);

sdfFModStM = cell(nU,nS); sdfFModStNMM = cell(nU,nS); % full model
sdfFModMM = cell(nU,nS); sdfFModNMM = cell(nU,nS);

for c = 1:size(sdfStM,1) % unit IDs
    for s = 1:size(sdfStM,2) % stimulus types 
        [~,sdfStM{c,s}] = io.triggeredAverage(trigIdxS{s},V(c,:),w,fS,0);
        [~,sdfStNMM{c,s}] = io.triggeredAverage(trigIdxSNM{s},V(c,:),w,fS,0);
        [~,sdfMM{c,s}] = io.triggeredAverage(trigIdxM{s},V(c,:),w,fS,0);
        [~,sdfNMM{c,s},tSDF] = io.triggeredAverage(trigIdxNM{s},V(c,:),w,fS,0);

        [~,sdfFModStM{c,s}] = io.triggeredAverage(trigIdxS{s},Vfl(c,:),w,fS,0);
        [~,sdfFModStNMM{c,s}] = io.triggeredAverage(trigIdxSNM{s},Vfl(c,:),w,fS,0);
        [~,sdfFModMM{c,s}] = io.triggeredAverage(trigIdxM{s},Vfl(c,:),w,fS,0);
        [~,sdfFModNMM{c,s}] = io.triggeredAverage(trigIdxNM{s},Vfl(c,:),w,fS,0);

    end
end
tSDF = tSDF*1000; % in ms

% save SDFs for neural data and model =====================================
if ~exist(fullfile(analDir,sesID),'dir')
    mkdir(fullfile(analDir,sesID))
end
cd(fullfile(analDir,sesID))
save('sdfData','sdfStM','sdfStNMM','sdfMM','sdfNMM','tSDF')
save('sdfFullModel','sdfFModStM','sdfFModStNMM','sdfFModMM','sdfFModNMM','tSDF')


% compute Motion Index ====================================================
rtVdt   = cell(nU,nS);
rtVdtNM   = cell(nU,nS);
for c = 1:nU % # units
    for s = 1:nS % # different stimulus types
        iSt = StM{s}(~isnan(StM{s}) & ~isnan(StNMM{s}));        
        for n = 1:length(iSt)
            its = ts>iSt(n)+tWinsc(1) & ts<=iSt(n)+tWinsc(2);
            rtVdt{c,s}(n)   = double(mean(V(c,its)));
        end

        iSt = StNMM{s}(~isnan(StM{s}) & ~isnan(StNMM{s}));
        for n = 1:length(iSt)
            its = ts>iSt(n)+tWinsc(1) & ts<=iSt(n)+tWinsc(2);
            rtVdtNM{c,s}(n)   = double(mean(V(c,its)));
        end
    end
end
%% average weighted by number of stimuli per condition (we had matched the
% stimulus values between movement conditions)
wMI_Vdt = nan(1,nU);
for c = 1:nU % units
    mOn = nanmean(cell2mat(rtVdt(c,:)));
    mOff = nanmean(cell2mat(rtVdtNM(c,:)));
    wMI_Vdt(c) = (mOn-mOff)/(mOn+mOff);
end

% save MIs ============================================================
cd(fullfile(analDir,sesID))
save('MIs','wMI_Vdt')

   
