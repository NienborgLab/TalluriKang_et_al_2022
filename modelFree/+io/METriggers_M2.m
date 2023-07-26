function METriggers_M2(sesDir,rootDir,outDir,wi)

%% input:
% sesDir            session directory
% rootDir           data root dir
% tsDir             rootDir of time-stamps
% outDir            output 
% wi                width of smoothing kernel


%% LOAD ME DATA 
sesDir_raw = strrep(sesDir,'M2_','ma_');
sesDir_raw = io.convertSessionID(sesDir_raw,rootDir);


% check availability of required files
cd(fullfile(rootDir,sesDir_raw))
if ~exist('faceME.mat','file') || ~exist('bodyME.mat','file')
    fprintf('%s %s ME not available\n',sesDir_raw)
    return
end
if ~exist('vid.mat','file')
    fprintf('%s vid not available\n',sesDir_raw)
    return
end

% load ME 
load('bodyME.mat','ME')
bME = ME;
load('faceME.mat','ME')
fME = ME;  clear ME
load('vid.mat','vid') % time-stamps for video data
tV = vid.t; clear vid
load('fNm.mat','fNm') %file names to assign hemispheres

% smooth ME
fME.av = conv(fME.av,ones(1,wi),'same')/wi;
bME.av = conv(bME.av,ones(1,wi),'same')/wi;

% equal lengths for both views 
vlen = [length(bME.av), length(fME.av)];
bME.av = bME.av(1:min(vlen));
fME.av = fME.av(1:min(vlen));

% exInfo
cd ..
ex = io.loadExInfo(sesDir_raw);
if isempty(ex)
    fprintf('%s exfile not available\n',sesDir_raw)
    return
end

% assign units to hemispheres
pat = '(?<probe>\w+).elec(?<electrode>\w+)\s+spk_(?<cluster>\w+)';
units = cellfun(@(x) regexp(x, pat, 'names'), fNm, 'uni', 1);

idx = find([ex.setup.gv.probe.leftHemisphereRecorded]==0);

if isfield(ex.setup.gv.probe(idx),'port')
    rPort   = ex.setup.gv.probe(idx).port;  % name of port for left visual field
else
    disp('no port recorded')
    rPort = 'A1';
end
hIdx{1} = find(ismember({units(:).probe},rPort)); 
rIdx    = find(~ismember({units(:).probe},rPort));
  
if ~isempty(rIdx)
    hIdx{2} = rIdx;  % indices for units with RF in left (1) right (2) visual field
end

%%
fprintf('processing... %s \n',sesDir)
% only use stimulus presentations on rewarded trials
tr = ex.Trials;
full_trials = find(abs([tr.Reward])>0 & [tr.instructionTrial]==0);

isi = []; % inter-stimulus interval (sec)
for n = 1:length(full_trials)
    itr = full_trials(n); 
    % get Isis
    if itr>1 && ~isempty(tr(itr-1).Start)
        Isi = tr(itr).Start(1)-tr(itr-1).Start(end);
    elseif itr>2 &&  ~isempty(tr(itr-2).Start)
        Isi = tr(itr(1)).Start(1)-tr(itr(1)-2).Start(end);
    else
        Isi = .5; % if we have 2 failed fixations assume IsI is >=.5
    end
    isi    = [isi,Isi]; %#ok<*AGROW>
end
tr = tr(full_trials);

% stimulus onset times (absolute time in session)
sOn = nan(length(tr),1);
for n = 1:length(tr)
    sOn(n) = tr(n).Start(1)-tr(n).TrialStart + tr(n).nevTrialStart;
end

%% match indices for stimulus onset times to video timestamp epochs
% use times from 0.2sec before to 2 sec after stimulus onset
idxC = cell(1,length(sOn));
for n = 1:length(tr)
    idxC{n} = find(tV>=sOn(n)-0.2 & tV<=sOn(n)+2);
end
idx = cell2mat(idxC);

sPer = zeros(1,length(bME.av));
sPer(idx) = true;

% identify high and low-motion periods based on percentiles
% [bqs] = quantile(bME.av(sPer==1),[.5, .8]);  % values used for M1 
% [fqs] = quantile(fME.av(sPer==1),[.5, .8]);  % values used for M1 
[bqs] = quantile(bME.av(sPer==1),[.2, .9]);  
[fqs] = quantile(fME.av(sPer==1),[.2, .9]);  
frameRate   = 20;
tol         = 0.5;
noMovTolf   = 0.2;
noMovTolb   = 0.3;
latTol      = 0.05; % match stimulus onset within 50ms

% face vid
idxMtf       = io.findMotion_simple(fME.av,frameRate,fqs(2),tol);
idxNMt_rawf  = io.findNoMotion_simple(fME.av,frameRate,fqs(1),noMovTolf);

% body vid
idxMtb       = io.findMotion_simple(bME.av,frameRate,bqs(2),tol);
idxNMt_rawb  = io.findNoMotion_simple(bME.av,frameRate,bqs(1),noMovTolb);

% use all possible no motion frames as no-motion onset
idxNMt = [];
for n = 1:size(idxNMt_rawf,1)
    idxNMt = [idxNMt, idxNMt_rawf(n,1):idxNMt_rawf(n,2)];
end
for n = 1:size(idxNMt_rawb,1)
    idxNMt = [idxNMt, idxNMt_rawb(n,1):idxNMt_rawb(n,2)];
end
idxNMt = unique(idxNMt)';

% use all possible motion frames, combined across views
idxM = [];
for n = 1:size(idxMtf,1)
    idxM = [idxM,idxMtf(n,1)-round(frameRate*noMovTolf):idxMtf(n,2)+round(frameRate*noMovTolf)];
end
for n = 1:size(idxMtb,1)
    idxM = [idxM,idxMtb(n,1)-round(frameRate*noMovTolb):idxMtb(n,2)+round(frameRate*noMovTolb)];
end
idxM = unique(idxM);
idx = ~ismember(idxNMt,idxM);
clear idxM idxNMt_raw
idxMt  = unique([idxMtf(:,1);idxMtb(:,1)]);
idxNMt = idxNMt(idx);

%
% get movement onset relative to stimulus onset, separated by stimulus and
% attention type
sCue = round([tr.Dc].*[tr.hdx]*100)/100;
sNCue = round(([tr.Dc2].*[tr.hdx2]+2)*100)/100; % add 2 to label uncued stimuli
att = sign([tr.x0]); % >0 attention on the Right, <0 att on the Left

sR = nan(1,length(sCue));  % right stimulus
sL = nan(1,length(sCue));  % left stimulus

sR(att>0) = sCue(att>0);
sR(att<0) = sNCue(att<0);

sL(att<0) = sCue(att<0);
sL(att>0) = sNCue(att>0);

% stimulus values, separated by hemifield
S{1} = sL;  % 1: left
if length(hIdx)>1 % we have units with RFs in both hemifields
    S{2} = sR;  % 2: right
end

tID = 1:length(full_trials); %unique ID for the included trials

%%
tMt     = tV(idxMt(:,1));
tNoMt   = tV(idxNMt(:,1));
cosStr{1} = '<=1';
cosStr{2} = '>1';

% initialize
ev = cell(length(hIdx),length(cosStr)); % align on all stimuli, separately for L/R stimulus
evM = cell(length(hIdx),length(cosStr));

dS = cell(length(hIdx),length(cosStr)); % lag between stimulus and movement onset
St = cell(length(hIdx),length(cosStr));  % stimulus onset time
StM = cell(length(hIdx),length(cosStr)); % stimulus onset time, matched for movement/no movement
tIDs = cell(length(hIdx),length(cosStr));  % trial ID
tIDsM = cell(length(hIdx),length(cosStr));

isiM = cell(length(hIdx),length(cosStr));  % ISI
isiNM = cell(length(hIdx),length(cosStr));
isiMM = cell(length(hIdx),length(cosStr));
isiNMM = cell(length(hIdx),length(cosStr));

SvalM = cell(length(hIdx),length(cosStr));  % stimulus value
SvalNM = cell(length(hIdx),length(cosStr));
SvalMM = cell(length(hIdx),length(cosStr));
SvalNMM = cell(length(hIdx),length(cosStr));

evNM = cell(length(hIdx),length(cosStr)); % event
evNMM = cell(length(hIdx),length(cosStr));

dSNM = cell(length(hIdx),length(cosStr)); % lag between stimulus and no-movement onset
StNM = cell(length(hIdx),length(cosStr));
StNMM = cell(length(hIdx),length(cosStr));

tIDsNM = cell(length(hIdx),length(cosStr)); % trial IDs
tIDsNMM = cell(length(hIdx),length(cosStr)); % trial IDs

IdxMt = cell(length(hIdx),length(cosStr));  % included ME-on indices
IdxNMt = cell(length(hIdx),length(cosStr)); % included ME-off indices
IdxMtM = [];  % included ME-on indices, matched stimuli
IdxNMtM = []; % included ME-off indices, matched stimuli

for h = 1:length(hIdx)  % hemifields
    for c = 1:length(cosStr)
        % select stimuli of specific contrast
        ix    = eval(['find(S{h}' cosStr{c} ');']);
        ix = ix(ix>1);
        istim = sOn(ix);
        % motion epochs 
        for n = 1:size(tMt,2)
            idx = find(istim>tMt(n)-2 & istim<tMt(n)+.5);  % select stimuli that occurred within 0.5s of a movement onset
            if ~isempty(idx)
                ev{h,c} = [ev{h,c},ones(1,length(idx))*tMt(n)]; % movement events
                dS{h,c} = [dS{h,c}, istim(idx)'-tMt(n)]; % lag beween stimulus and movement
                St{h,c} = [St{h,c}, istim(idx)']; % stimulus times
                isiM{h,c} = [isiM{h,c},isi(ix(idx))]; % isi of current stimuli
                SvalM{h,c} = [SvalM{h,c},S{h}(ix(idx))];
                tIDs{h,c} = [tIDs{h,c},tID(ix(idx))];
                for n2 = 1:length(idx)% stimulus aligned motion
                    [~,tIdx] = min(abs(tV-istim(idx(n2)))); 
                    IdxMt{h,c} = [IdxMt{h,c}, tIdx];
                end
            end      
        end
        % no motion epochs and window +/- 0.5sec around it
        for n = 1:size(tNoMt,2)
            idx = find(istim>tNoMt(n)-2 & istim<tNoMt(n)+0.5);
            if ~isempty(idx)
                evNM{h,c} = [evNM{h,c},ones(1,length(idx))*tNoMt(n)];
                dSNM{h,c} = [dSNM{h,c}, istim(idx)'-tNoMt(n)];
                StNM{h,c} = [StNM{h,c}, istim(idx)'];
                isiNM{h,c} = [isiNM{h,c},isi(ix(idx))];
                SvalNM{h,c} = [SvalNM{h,c},S{h}(ix(idx))];
                tIDsNM{h,c} = [tIDsNM{h,c},tID(ix(idx))];
                
                for n2 = 1:length(idx)
                    [~,tIdx] = min(abs(tV-istim(idx(n2))));
                    IdxNMt{h,c} = [IdxNMt{h,c}, tIdx];
                end

            end
        end

        % match no motion onset times to motion onset
        for n = 1:length(dS{h,c})
            iDS = dS{h,c}(n);
            iisi = isiM{h,c}(n);
            iSval = SvalM{h,c}(n);
            % match stimulus, stimulus onset within 50ms, ISI within 300ms
             idx = find(abs(dSNM{h,c}-iDS)<latTol & abs(isiNM{h,c}-iisi)<.3 & SvalNM{h,c}==iSval); 

            if ~isempty(idx)
                evNMM{h,c} = [evNMM{h,c},evNM{h,c}(idx)];
                evM{h,c}   = [evM{h,c},ev{h,c}(n)*ones(1,length(idx))];
                StM{h,c}   = [StM{h,c},St{h,c}(n)*ones(1,length(idx))];
                StNMM{h,c} = [StNMM{h,c},StNM{h,c}(idx)];
                isiMM{h,c} = [isiMM{h,c},isiM{h,c}(n)*ones(1,length(idx))];
                isiNMM{h,c} = [isiNMM{h,c},isiNM{h,c}(idx)];
                IdxMtM   = [IdxMtM,IdxMt{h,c}(n)*ones(1,length(idx))];
                IdxNMtM    = [IdxNMtM,IdxNMt{h,c}(idx)];
                SvalMM{h,c} = [SvalMM{h,c},SvalM{h,c}(n)*ones(1,length(idx))];
                SvalNMM{h,c} = [SvalNMM{h,c},SvalNM{h,c}(idx)];
                
                tIDsM{h,c} = [tIDsM{h,c},tIDs{h,c}(n)*ones(1,length(idx))];
                tIDsNMM{h,c} = [tIDsNMM{h,c},tIDsNM{h,c}(idx)];
             else % remove epoch if we cannot match up the stimulus
                    %fprintf('not possible to match %d %d\n',c,n)
                    dS{h,c}(n) = NaN;
                    ev{h,c}(n) = nan;
                    St{h,c}(n) = nan;
                    IdxMt{h,c}(n) = nan;
                    
            end

        end

    end
end

%% remove stimuli included in both sets and match by stimulus and attention 
%% condition again
for h = 1:size(StM,1)
    for c = 1:size(StM,2)
        [~,idx] = setdiff(StM{h,c},StNMM{h,c});
        StM{h,c} = StM{h,c}(idx);
        SvalMM{h,c} = SvalMM{h,c}(idx);
        tIDsM{h,c} = tIDsM{h,c}(idx);

        [~,idx] = setdiff(StNMM{h,c},StM{h,c});
        StNMM{h,c} = StNMM{h,c}(idx);
        SvalNMM{h,c} = SvalNMM{h,c}(idx);
        tIDsNMM{h,c} = tIDsNMM{h,c}(idx);
    end
end

% ensure we have the same number of stimuli in each set, and we match the
% type of stimulus for each included epoch
for h = 1:size(StM,1)
    for c = 1:size(StM,2)
        istM = [];
        istNM = [];
        iSvalM = [];
        iSvalNM = [];
        
        itIDsM = []; % trialIDs
        itIDsNM = [];

        for s = 1:length(StM{h,c})
            iSval = SvalMM{h,c}(s); % current stimulus value
            idx = find(SvalNMM{h,c}==iSval);
            
            if ~isempty(idx)
                istM = [istM,StM{h,c}(s)]; 
                iSvalM = [iSvalM,iSval];
                itIDsM = [itIDsM,tIDsM{h,c}(s)];
                    
                [~,iidx] = setdiff(StNMM{h,c}(idx),istNM);
                if ~isempty(iidx)
                    rp = randperm(length(iidx));
                    istNM = [istNM,StNMM{h,c}(idx(iidx(rp(1))))];
                    iSvalNM = [iSvalNM,SvalNMM{h,c}(idx(iidx(rp(1))))];
                    itIDsNM  = [itIDsNM,tIDsNMM{h,c}(idx(iidx(rp(1))))];
                else
                    rp = randperm(length(idx));
                    istNM = [istNM,StNMM{h,c}(idx(rp(1)))];
                    iSvalNM = [iSvalNM,SvalNMM{h,c}(idx(rp(1)))];
                    itIDsNM = [itIDsNM,tIDsNMM{h,c}(idx(rp(1)))];
                end
            end
        end
        StM{h,c} = istM;
        StNMM{h,c} = istNM;
        SvalMM{h,c} = iSvalM;
        SvalNMM{h,c} = iSvalNM;
        
        tIDsM{h,c} = itIDsM;
        tIDsNMM{h,c} = itIDsNM;
    end
end

% save trigger times and preprocessed data ================================
if ~exist(fullfile(outDir,sesDir),'dir')
    mkdir(fullfile(outDir,sesDir))

end
cd(fullfile(outDir,sesDir))

save('both_trigTimes','StM','StNMM','evM','evNMM',...
     'tIDsM','tIDsNMM','SvalMM','SvalNMM','sOn','S') 

save('unitByHemisphere','hIdx')  % keep track of hemisphere information
