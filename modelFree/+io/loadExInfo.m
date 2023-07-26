function ex = loadExInfo(fDir)
% function ex = loadExInfo(fDir)
% convenience function to load exInfo in given dir

curDir                  = cd(fDir);
exIFn                   = dir('*exInfo.mat');
load(exIFn(1).name);

cd(curDir)