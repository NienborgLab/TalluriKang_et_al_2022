function METriggers_M1(sesDir,rootDir,outDir,wi)

%% input:
% sesDir            session directory
% rootDir           data root dir
% tsDir             rootDir of time-stamps
% outDir            output 
% wi                width of smoothing kernel


%% LOAD ME DATA 
sesDir_raw = strrep(sesDir,'M1_','le_');
sesDir_raw = io.convertSessionID(sesDir_raw,rootDir);

% check availability of required files
cd(fullfile(rootDir,sesDir_raw))
if ~exist('faceME.mat','file') || ~exist('bodyME.mat','file')
    fprintf('%s %s', 'ME not available\n',sesDir_raw)
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

%%
fprintf('processing... %s \n',sesDir)

% only use stimulus presentations on rewarded trials, i.e. for which
% we have 4 consecutive valid stimulus presentations
tr = ex.Trials;
full_trials = find([tr.RewardSize]>0);
full_trials = full_trials(full_trials>3);

idx = []; % trial index
rnk = []; % stimulus rank
isi = []; % inter-stimulus interval (sec)
for n = 1:length(full_trials)
    itr     = full_trials(n)-3:full_trials(n);
    if any([tr(itr).Reward]==0)
        continue
    end
    Isi     = nan(1,3);
    for s = 1:3
        Isi(s) = tr(itr(s+1)).Start(1)-tr(itr(s)).Start(end);
    end
    idx     = [idx, itr]; %#ok<*AGROW>
    rnk     = [rnk, 1:4]; %#ok<*AGROW>
    
    % get Isis to preceding completed trial, 
    % check for up to 2 preceding stimulus presentations, default is 0.5sec 
    if itr(1)>1 && ~isempty(tr(itr(1)-1).Start)
        Isi = [tr(itr(1)).Start(1)-tr(itr(1)-1).Start(end),Isi];
    elseif itr(1)>2 &&  ~isempty(tr(itr(1)-2).Start)
        Isi = [tr(itr(1)).Start(1)-tr(itr(1)-2).Start(end),Isi];
    else
        Isi = [.5 Isi]; %#ok<*AGROW>
    end
    isi     = [isi,Isi]; %#ok<*AGROW>
end
tr = tr(idx);

% stimulus onset times (absolute time in session)
sOn = nan(length(tr),1);
for n = 1:length(tr)
    sOn(n) = tr(n).Start(1)-tr(n).TrialStart + tr(n).nevTrialStart;
end

%% match indices for stimulus onset times to video timestamp epochs
idxC = cell(1,length(sOn));
for n = 1:length(tr)
    % 200ms before to 500ms after stimulus onset
    idxC{n} = find(tV>=sOn(n)-0.2 & tV<=sOn(n)+0.5);
end
idx = cell2mat(idxC);

sPer = zeros(1,length(fME.av));
sPer(idx) = true;
%%
% identify high and low-motion periods based on percentiles during stimulus
% presentations
[fqs] = quantile(fME.av(sPer==1),[.5, .8]);
[bqs] = quantile(bME.av(sPer==1),[.5, .8]);
frameRate   = 20;
tol         = 0.5;   % merge movement epochs occuring within this time window (sec)
noMovTolf   = 0.2;   % require this amount of time without movement (sec); face view
noMovTolb   = 0.3;   % require this amount of time without movement (sec); body view

% face view
idxMtf       = io.findMotion_simple(fME.av,frameRate,fqs(2),tol);
idxNMt_rawf  = io.findNoMotion_simple(fME.av,frameRate,fqs(1),noMovTolf);

% body view
idxMtb       = io.findMotion_simple(bME.av,frameRate,bqs(2),tol);
idxNMt_rawb  = io.findNoMotion_simple(bME.av,frameRate,bqs(1),noMovTolb);

% use all possible no motion frames as no-movement onset, i.e. take all
% frames between no-movement onset and no-movement offset
idxNMt = [];
for n = 1:size(idxNMt_rawf,1)
    idxNMt = [idxNMt, idxNMt_rawf(n,1):idxNMt_rawf(n,2)];
end
for n = 1:size(idxNMt_rawb,1)
    idxNMt = [idxNMt, idxNMt_rawb(n,1):idxNMt_rawb(n,2)];
end
idxNMt = unique(idxNMt)';

% use all possible motion frames, combined across views to identify
% motion frames, and exclude them from no-movement sets
idxM = [];
for n = 1:size(idxMtf,1)
    idxM = [idxM,idxMtf(n,1)-round(frameRate*noMovTolf):idxMtf(n,2)+...
        round(frameRate*noMovTolf)];
end
for n = 1:size(idxMtb,1)
    idxM = [idxM,idxMtb(n,1)-round(frameRate*noMovTolb):idxMtb(n,2)+...
        round(frameRate*noMovTolb)];
end
idxM = unique(idxM);
idx = ~ismember(idxNMt,idxM); % all no-movement onsets not included in movement onsets

idxMt  = unique([idxMtf(:,1);idxMtb(:,1)]);  %unique motion onset indices
idxNMt = idxNMt(idx);                        % possible no-motion onset indices

clear idxM idxNMt_raw bME fME idxNMt_rawb idxNMt_rawf idxMtb idxMtf idx

%%
% get movement onset relative to stimulus onset, separated by contrast level
[tr([tr.co]>100).co] = deal(0);  % set blank stimulus to 0% contrast
Cos = unique([tr.co]);
Cos = sort(Cos);

tMt     = tV(idxMt(:,1)); %times of movement onset (sec)
tNoMt   = tV(idxNMt(:,1)); %times of movement offset (sec)

ev = cell(1,length(Cos)); % movement events
evM = cell(1,length(Cos)); % movement events (matched)
evM_u = cell(1,length(Cos)); % movement events (matched), unique stimuli

dS = cell(1,length(Cos)); % lag between stimulus and motion onset
St = cell(1,length(Cos)); % stimulus times

StM = cell(1,length(Cos)); % stimulus times with movement

copM = cell(1,length(Cos)); % contrast of preceding stimulus w movement
copNM = cell(1,length(Cos)); % contrast of preceding stimulus w/o movement

% matched conditions for movement and no movement epochs
copMM = cell(1,length(Cos)); % contrast of preceding stimulus w movement (matched)
copNMM = cell(1,length(Cos)); % contrast of preceding stimulus w/o movement (matched)

rnkM = cell(1,length(Cos)); % rank of stimulus, w movement
rnkNM = cell(1,length(Cos));  % rank of stimulus, w/o movement
rnkMM = cell(1,length(Cos));  % rank of stimulus, w movement (matched)
rnkNMM = cell(1,length(Cos)); % rank of stimulus, w/o movement (matched)

isiM = cell(1,length(Cos));  % isi for movement epochs
isiNM = cell(1,length(Cos)); % isi for no movement epochs
isiMM = cell(1,length(Cos)); % isi for movement epochs (matched)
isiNMM = cell(1,length(Cos)); % isi for no movement epochs (matched)

evNM = cell(1,length(Cos)); % no movement events
evNMM = cell(1,length(Cos)); % no movement events (matched)
evNMM_u = cell(1,length(Cos)); % no movement events (matched), unique stimuli

dSNM = cell(1,length(Cos)); % lag between stimulus and no movement onset
StNM = cell(1,length(Cos)); % stimuli w/o movement
StNMM = cell(1,length(Cos)); % stimuli w/o movement (matched)
IdxMt = cell(1,length(Cos));  % included movement-on indices
IdxNMt = cell(1,length(Cos)); % included no movement-on indices
IdxMtM = [];  % included movement-on indices, matched stimuli
IdxNMtM = []; % included no-movement-on indices, matched stimuli

for c = 1:length(Cos)
    % select stimuli of specific contrast
    ix    = find([tr.co]==Cos(c));
    ix = ix(find(ix>1));
    istim = sOn(ix);  % stimuli of one type
    
    % identify motion epochs and stimuli in window +/-0.7sec around it
    for n = 1:size(tMt,2)
        idx = find(istim>tMt(n)-.7 & istim<tMt(n)+.7);  % select stimuli 
        %                         that occured within 0.7s of a movement
        if ~isempty(idx)
            ev{c} = [ev{c},ones(1,length(idx))*tMt(n)]; % movement events
            dS{c} = [dS{c}, istim(idx)'-tMt(n)]; % lag beween stimulus and movement
            St{c} = [St{c}, istim(idx)']; % stimulus times
            copM{c} = [copM{c}, [tr(ix(idx)-1).co]]; % contrast of perceding stimulus
            rnkM{c} = [rnkM{c},rnk(ix(idx))]; % rank of current stimlus
            isiM{c} = [isiM{c},isi(ix(idx))]; % isi of current stimuli
            
            for n2 = 1:length(idx)% ME index for stimuli within movement window
                [~,tIdx] = min(abs(tV-istim(idx(n2)))); 
                IdxMt{c} = [IdxMt{c}, tIdx];
            end
        end      
    end
    
    % no motion epochs and stimuli in window +/-0.5sec around it
    for n = 1:size(tNoMt,2)
        idx = find(istim>tNoMt(n)-0.5 & istim<tNoMt(n)+0.5);
        if ~isempty(idx)
            evNM{c} = [evNM{c},ones(1,length(idx))*tNoMt(n)];
            dSNM{c} = [dSNM{c}, istim(idx)'-tNoMt(n)];
            StNM{c} = [StNM{c}, istim(idx)'];
            rnkNM{c} = [rnkNM{c},rnk(ix(idx))];
            isiNM{c} = [isiNM{c},isi(ix(idx))];
            copNM{c} = [copNM{c},[tr(ix(idx)-1).co]];
            for n2 = 1:length(idx)
                [~,tIdx] = min(abs(tV-istim(idx(n2))));
                IdxNMt{c} = [IdxNMt{c}, tIdx];
            end

        end
    end
    
    % now match conditions for no motion onset epochs to motion onset
    % events
    for n = 1:length(dS{c})
        iDS = dS{c}(n);
        ico = copM{c}(n);
        irnk = rnkM{c}(n);
        iisi = isiM{c}(n);
        
        % match stimulus onset within 10ms, isi within 50ms, rank and
        % preceding stimulus
         idx = find(abs(dSNM{c}-iDS)<.01 & abs(isiNM{c}-iisi)<.05 & ...
             copNM{c}==ico & rnkNM{c}==irnk); 

        if ~isempty(idx)
            evNMM{c} = [evNMM{c},evNM{c}(idx)];
            evM{c}   = [evM{c},ev{c}(n)*ones(1,length(idx))]; % duplicate movement events to match
            StM{c}   = [StM{c},St{c}(n)*ones(1,length(idx))]; % duplicate stimuli to match
            StNMM{c} = [StNMM{c},StNM{c}(idx)]; % include only matched stimuli
            copNMM{c} = [copNMM{c},copNM{c}(idx)]; %
            copMM{c} = [copMM{c},ico*ones(1,length(idx))];    
            
            rnkNMM{c} = [rnkNMM{c},rnkNM{c}(idx)];
            rnkMM{c} = [rnkMM{c},irnk*ones(1,length(idx))];
            isiMM{c} = [isiMM{c},isiM{c}(n)*ones(1,length(idx))];
            isiNMM{c} = [isiNMM{c},isiNM{c}(idx)];
            IdxMtM   = [IdxMtM,IdxMt{c}(n)*ones(1,length(idx))];
            IdxNMtM    = [IdxNMtM,IdxNMt{c}(idx)];
        end
    end
end
%%   
%% remove stimuli included in both sets, then match by rank, ISI 
%% and previous CO again
for c = 1:size(StM,2)
    [~,idx] = setdiff(StM{c},StNMM{c});
    StM{c}   = StM{c}(idx);
    rnkMM{c} = rnkMM{c}(idx);
    isiMM{c} = isiMM{c}(idx);
    copMM{c} = copMM{c}(idx);
    evM_u{c}  = evM{c}(idx);
    
    [~,idx] = setdiff(StNMM{c},StM{c});
    StNMM{c} = StNMM{c}(idx);
    rnkNMM{c} = rnkNMM{c}(idx);
    isiNMM{c} = isiNMM{c}(idx);
    copNMM{c} = copNMM{c}(idx);
    evNMM_u{c}  = evNMM{c}(idx);
end
for c = 1:size(StM,2)
    istM = [];
    istNM = [];
    irnkM = [];
    irnkNM = [];
    ievM = [];
    ievNM =[];
    for s = 1:length(StM{c})
        irn = rnkMM{c}(s);
        ico = copMM{c}(s);
        IsI = isiMM{c}(s);
        
        % match by rank, previous contrast, ISI (within 100ms)
        idx = find(rnkNMM{c}==irn & copNMM{c} ==ico & abs(isiNMM{c}-IsI)<0.1);
        if ~isempty(idx)
            istM = [istM,StM{c}(s)];
            irnkM = [irnkM,rnkMM{c}(s)];
            ievM = [ievM,evM_u{c}(s)];
            [~,iidx] = setdiff(StNMM{c}(idx),istNM);
            if ~isempty(iidx)
                rp = randperm(length(iidx));
                istNM = [istNM,StNMM{c}(idx(iidx(rp(1))))];
                irnkNM = [irnkNM,rnkNMM{c}(idx(iidx(rp(1))))];
                ievNM = [ievNM,evNMM_u{c}(idx(iidx(rp(1))))];
            else % randomly select one of the stimuli we already included
                rp = randperm(length(idx));
                istNM = [istNM,StNMM{c}(idx(rp(1)))];
                irnkNM = [irnkNM,rnkNMM{c}(idx(rp(1)))];
                ievNM = [ievNM,evNMM_u{c}(idx(rp(1)))];
            end
        end
    end
    % only use stimuli and epochs that could be matched
    StM{c} = istM;  
    StNMM{c} = istNM;  % no movement stimuli matched in rank, preceding CO and ISI
    rnkMM{c} = irnkM;
    rnkNMM{c} = irnkNM;
    evM_u{c} = ievM; % only one event per stimulus included
    evNMM_u{c} = ievNM; % only one event per stimulus included
end

% save trigger times and preprocessed data ================================
if ~exist(fullfile(outDir,sesDir),'dir')
    mkdir(fullfile(outDir,sesDir))

end
cd(fullfile(outDir,sesDir))
save('both_trigTimes','StM','StNMM','evM','evNMM','rnkMM','rnkNMM')



