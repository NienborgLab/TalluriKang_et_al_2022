ex_names = {'M1_0008', 'M1_0010', 'M1_0012', 'M1_0013', 'M1_0014', 'M1_0031',...
    'M1_0032', 'M1_0033', 'M1_0034', 'M1_0036', 'M1_0037', 'M1_0038', 'M1_0039',...
    'M1_0040', 'M1_0042', 'M1_0043', 'M1_0044', 'M1_0045', 'M1_0046', 'M1_0048',...
    'M1_0049', 'M1_0050_all1', 'M1_0050_all2', 'M1_0052', 'M1_0054', 'M1_0055',...
    'M1_0056', 'M1_0057', 'M1_0058_15', 'M1_0058_all', 'M1_0059', 'M1_0060',...
    'M1_0063', 'M1_0064', 'M1_0065', 'M1_0066_11', 'M1_0066_all', 'M1_0067_all1',...
    'M1_0067_all2', 'M1_0068', 'M1_0069' 'M1_0071', 'M1_0072_all1', 'M1_0072_all2',...
    'M1_0073_all1', 'M1_0073_all2', 'M1_0074', 'M1_0075', 'M1_0076', 'M1_0077',...
    'M1_0078_all1', 'M1_0078_all2', 'M1_0079', 'M1_0080', 'M2_0803', 'M2_0810',...
    'M2_0822', 'M2_0823', 'M2_0824'};

addpath preprocess/
addpath analysis/
addpath tools/
addpath visualise/

clc;close all;rng(1);
myfigureprops;
% we will use a 300 ms window for every trial before stimulus onset
time_pre_stim = 0.3; % in seconds
% we will use a bin-width of 16.67 ms, matching the stimulus display
% frequency
dt = 1/60;
% we can also specify the maximum trial duration using the variable
% time_per_trial: for this analysis, we will a time-window that includes
% the stimulus presentation window (2 s), and a max 1 s post-stimulus
% window
%% get some data from example recording for sample preprocessing
raw_datadir = '../data/raw/';
ex_name = 'M1_0060'; % another example session with raw data is 'M2_0822'
if strcmp(ex_name(1:2), 'M1')
    [D, D_eye] = M1.preprocess_dataset(raw_datadir, ex_name , 'time_pre_stim', time_pre_stim, 'dt', dt);
elseif strcmp(ex_name(1:2), 'M2')
    [D, D_eye] = M2.preprocess_dataset(raw_datadir, ex_name , 'time_pre_stim', time_pre_stim, 'dt', dt);
end
save(sprintf('../data/preprocessed/%s/preprocessed_data.mat', ex_name), 'D', 'D_eye', '-v7.3');
clear D D_eye

%% main analysis

for ex = 1:length(ex_names)
    fprintf('%s \n', ex_names{ex});
    exname = ex_names{ex};
    linear_Model(exname);
end
harvest_vars();
fig2('all');
fig2('M1');
fig2('M2');
figSI();

% run models with eye
for ex = 1:length(ex_names)
    fprintf('%s \n', ex_names{ex});
    exname = ex_names{ex};
    linear_Model(exname, 'include_eye', 1);
end
harvest_vars('include_eye', 1);
fig2('all', 'include_eye', 1);
figSI('include_eye', 1);