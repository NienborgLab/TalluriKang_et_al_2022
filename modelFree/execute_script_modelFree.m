addpath(pwd)
cd .. ; cd ..
rootDir     =    pwd; 
rawDataDir  =    fullfile(rootDir, '\data\raw');
trigDir     =    fullfile(rootDir, '\data\modelFree\preprocessed\triggerTimes');
algTrigDir  =    fullfile(rootDir, '\data\modelFree\preprocessed\alignedTriggers');
analDir     =    fullfile(rootDir, '\data\modelFree\analysis');
tsDir       =    fullfile(rootDir, '\data\preprocessed');  
rrmDir      =    fullfile(rootDir, '\data\analysis'); 
figDir      =    fullfile(rootDir, '\data\results');

%% preprocessing:
% create list of sessions
cd(tsDir)
lst = dir('M*_*');  % list of all experimental sessions

%% preprocess the data to get the onsets (triggerTimes) of high vs low movement
for n = 1:length(lst)
    sesNm = lst(n).name; % session name
    if contains(sesNm,'M1_')
        wi = 12;
        io.METriggers_M1(sesNm,rawDataDir,trigDir,wi);
    elseif contains(sesNm,'M2_')
        wi = 10;
        io.METriggers_M2(sesNm,rawDataDir,trigDir,wi);
    end
end

%% align trigger-times to time stamps of preprocessed data for rrm and fits
for n = 1:length(lst)  
    io.alignTriggers(lst(n).name,trigDir,tsDir,algTrigDir);
end


%% analyzing: compute the triggered PSTHs, movement index, attention index
for n = 1:length(lst) 
    sesNm = lst(n).name; % session name
    if contains(sesNm,'M1_')
        io.trigSDFandMI_M1(sesNm,tsDir,trigDir,algTrigDir,rrmDir,analDir); 
    elseif contains(sesNm,'M2_')
        io.trigSDFandMI_M2(sesNm,tsDir,trigDir,algTrigDir,rrmDir,analDir); 
    end
end

% harvest and plot the data
plot.fig3A(analDir,rrmDir,tsDir,figDir)
plot.fig3B(analDir,rrmDir,tsDir,figDir)
plot.supplFigS6(analDir,rrmDir,tsDir,figDir)
