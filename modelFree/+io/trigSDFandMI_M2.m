function trigSDFandMI_M2(sesID,tsDir,trigDir,algTrigDir,rrmDir,analDir) 

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

fS = 60; % sampling frequency of data/model (60Hz)
tWin = [-1 2]; % timewindow for SDF
w = abs(tWin*fS); % width window around trigger of SDF
tWinsc = [0.15 2]; % timewindow for spike-counts to compute MotionIndex (MI)
% and attention index (AI); exclude onset transient

tic

fprintf('processing... %s \n',sesID)
% load triggers
load(fullfile(algTrigDir,sesID,'both_trigIdx'),...
    'trigIdxS','trigIdxSNM','trigIdxM','trigIdxNM');

% load hemisphere info
load(fullfile(trigDir,sesID,'unitByHemisphere'),'hIdx');

% load trigger times
load(fullfile(trigDir,sesID,'both_TrigTimes'),'StM','StNMM','SvalMM','SvalNMM','sOn','S');
% sOn: all stimulus onset times 
% S: all stimulus values in L/R hemifield

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
sdfStM = cell(nU,nS); sdfStNMM = cell(nU,nS);
sdfMM = cell(nU,nS); sdfNMM = cell(nU,nS);

sdfFModStM = cell(nU,nS); sdfFModStNMM = cell(nU,nS);
sdfFModMM = cell(nU,nS); sdfFModNMM = cell(nU,nS);

for h = 1:length(hIdx)
    idx = hIdx{h};
    for c = 1:idx % unit IDs
        for s = 1:size(sdfStM,2) % stimulus types 
            [~,sdfStM{idx(c),s}] = io.triggeredAverage(trigIdxS{h,s},V(idx(c),:),w,fS,0);
            [~,sdfStNMM{idx(c),s}] = io.triggeredAverage(trigIdxSNM{h,s},V(idx(c),:),w,fS,0);
            [~,sdfMM{idx(c),s}] = io.triggeredAverage(trigIdxM{h,s},V(idx(c),:),w,fS,0);
            [~,sdfNMM{idx(c),s},tSDF] = io.triggeredAverage(trigIdxNM{h,s},V(idx(c),:),w,fS,0);

            [~,sdfFModStM{idx(c),s}] = io.triggeredAverage(trigIdxS{h,s},Vfl(idx(c),:),w,fS,0);
            [~,sdfFModStNMM{idx(c),s}] = io.triggeredAverage(trigIdxSNM{h,s},Vfl(idx(c),:),w,fS,0);
            [~,sdfFModMM{idx(c),s}] = io.triggeredAverage(trigIdxM{h,s},Vfl(idx(c),:),w,fS,0);
            [~,sdfFModNMM{idx(c),s}] = io.triggeredAverage(trigIdxNM{h,s},Vfl(idx(c),:),w,fS,0);

        end
    end
end
tSDF = tSDF*1000; % in ms

% save SDFs for neural data and model =====================================
dirNm = strrep(sesID,'ma_','M2_');
if ~exist(fullfile(analDir,dirNm),'dir')
    mkdir(fullfile(analDir,dirNm))
end
cd(fullfile(analDir,dirNm))
save('sdfData','sdfStM','sdfStNMM','sdfMM','sdfNMM','tSDF')
save('sdfFullModel','sdfFModStM','sdfFModStNMM','sdfFModMM',...
    'sdfFModNMM','tSDF')

% compute Motion Index and Attention Index ================================
rtVdt   = cell(nU,nS);
rtVdtNM   = cell(nU,nS);

% get unique stimulus values to compute attention index
% SvalMM: 2-by-2: hemisphere by attention condition
if size(SvalMM,1)>1
    uS = unique([cell2mat(SvalMM(1,:)), cell2mat(SvalMM(2,:))]);
else 
    uS = unique(cell2mat(SvalMM(1,:)));
end

for h = 1:length(hIdx) % hemisphere
    idx = hIdx{h};
    for c = 1:length(idx) % unit
        for s =1:length(uS) % stimulus
            
            % movement epochs
            iSt = cell2mat(StM(h,:));
            iStN = cell2mat(StNMM(h,:));
            ids = find(cell2mat(SvalMM(h,:)) == uS(s) & ~isnan(iSt) &~isnan(iStN));

            for n = 1:length(ids)
                its = ts>iSt(ids(n)) ...
                    +tWinsc(1) &  ts<=iSt(ids(n))+tWinsc(2);
                rtVdt{idx(c),s}(n) = double(mean(V(idx(c),its)));
            end

            % no movement epochs
            ids = find(cell2mat(SvalNMM(h,:)) == uS(s) & ~isnan(iSt) &~isnan(iStN));

            for n = 1:length(ids)
                its = ts>iStN(ids(n)) ...
                    +tWinsc(1) &  ts<=iStN(ids(n))+tWinsc(2);
                rtVdtNM{idx(c),s}(n) = double(mean(V(idx(c),its)));
            end
        end % stimulus
    end % unit
end % hemisphere

% average weighted by number of stimuli per condition (we match the
% stimulus values between conditions)
wMI_Vdt = nan(1,nU);
for c = 1:nU % units
    mOn = nanmean(cell2mat(rtVdt(c,:)));
    mOff = nanmean(cell2mat(rtVdtNM(c,:)));
    wMI_Vdt(c) = (mOn-mOff)/(mOn+mOff);
end

%% compute attention index ================================================
spRtVdt = cell(1,nU);
for c = 1:nU % units
    for n = 1:length(sOn)
        idx        = ts>sOn(n)+tWinsc(1) ...
            &  ts<=sOn(n)+tWinsc(2);
        spRtVdt{c}(n) = mean(V(c,idx));

    end
end

sCues   = uS(uS<1); % cued stimulus 
sNCues  = uS(uS>1); % uncued stimulus
AI_Vdt = nan(1,nU);
for h = 1:length(S) % hemisphere
    idx = hIdx{h};
    for c = 1:length(idx) % units
        aiByStimVdt = nan(1,length(sCues));
        for s = 1:length(sCues) % stimulus types
            aiVdt = mean(spRtVdt{idx(c)}(S{h}==sCues(s)));
            aoVdt = mean(spRtVdt{idx(c)}(S{h}==sNCues(s)));
            aiByStimVdt(s) = (aiVdt-aoVdt)/(aiVdt+aoVdt);
        end
        AI_Vdt(idx(c)) = nanmean(aiByStimVdt); 
    end
end

% save MIs and AIs ========================================================
cd(fullfile(analDir,sesID))
save('AI_MIs','AI_Vdt','wMI_Vdt')

   
