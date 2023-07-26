function [vSVD,SVD,ME] = motionSVD(videoFile, roiFile,nt,ns)

% This script is a modification of
% https://github.com/MouseLand/facemap/blob/main/matlab/computeSVDmotion.m
% to incooperate the data structure produced in the Nienborg lab.
%
% The three MATLAB scripts, svdecon.m, gather_try.m, normc.m, used in this
% script are % adopted from 
% https://github.com/MouseLand/facemap/tree/main/matlab/utils
% without modifications.
%
% Example:
% [vSVD,SVD] = motionSVD({'videoFile1.avi','videoFile2.avi'},...
%    {'roiFile.mat','roiFile.mat'},3,2);
%
% Inputs:
%
%   videoFile:  names of video files, cell(1,number of video files) 
%
%   roiFile:    names of ROI files, cell(1,number of video files), ROI
%               files are expected to be the same for a session.
%               A roiFile is a logical array of the same size of the video
%               frame in which pixels included in the analysis are set true.
%
%   nt:         temporal downsamling factor: video will be sampled at 1/nt
%               of the frame rate of the original video.
%
%   ns:         spatial downsamling factor: video size will be reduced to
%               1/ns^2 of the original size.
%
% Outputs:
%
%   vSVD.uME:   normalized left signular vectors of motion energy
%   vSVD.svdME: projections of signular vectors
%   vSVD.ROI:   reigion of interest in the video
%
%   SVD.uMot:   concantnated singular vectors
%   SVD.uME:    left singular vectors of SVD.uMot
%   SVD.sME:    singular values of SVD.uMot
%   SVD.vME:    right singular vectors of SVD.uMot
%
%   ME:         "motion energy" (abs of pixel difference)


if nargin == 2
    nt = 1;
    ns = 1;
elseif nargin == 3
    ns = 1;
end

useParallel = true;

nVideoFiles = length(videoFile);

vobj = cell(1,nVideoFiles);
nFrames = zeros(1,nVideoFiles);
for i=1:nVideoFiles
    vobj{i} = VideoReader(videoFile{i});

    nFrames(i) = vobj{i}.NumFrames;
    
    load(roiFile{i},'ROI')
    if i == 1
        roi = zeros(vobj{i}.Height,vobj{i}.Width,length(videoFile),'logical');
    end
    roi(:,:,i) = ROI;
end

if nVideoFiles > 1
    a = mean(single(roi),3);
    if ~isempty(setdiff(unique(a(:)),[0,1]'))
        error('video files might have different ROIs')
    end
    clear a
end
roi = roi(:,:,1);

% downsample roi if necessay
if ns > 1
    nh = floor(vobj{1}.Height/ns);
    nw = floor(vobj{1}.Width/ns);

    roi = single(roi);
    roi = squeeze(mean(mean(reshape(roi(1:nh*ns,1:nw*ns,1),...
        ns, nh, ns, nw), 1),3));

    roi(roi<0.5) = 0;
    roi = logical(roi);
end

% generate video indices for batch process
batchSize = 2000 * nt;
a = cell(nVideoFiles,1);
nBatches = zeros(nVideoFiles,1);
for i=1:nVideoFiles
    idxStart = (1:batchSize:nFrames(i))';
    idxEnd = [idxStart(2:end)-1;nFrames(i)];
    nBatches(i) = length(idxStart);
    a{i} = [ones(nBatches(i),1) * i,...
        idxStart, idxEnd,...
        idxEnd - idxStart + 1];
end
a = cell2mat(a);
indices = table;
indices.i = a(:,1);
indices.index = a(:,2:3);
indices.nFrames = a(:,4);
indices = table2struct(indices);

% poolSize needs to be adjucted depending on video size and length to avoid
% depletion of memory
if useParallel
    poolSize = 4;
    p = gcp('nocreate');
    if isempty(p)
        parpool(poolSize);
    else
        if p.NumWorkers ~= poolSize
            delete(p);
            parpool(poolSize);
        end
    end
end

% estimate average pixel values of ROI in the video and motion energy
nBatch = length(indices);
nPixels = sum(roi(:));
avgFrame = single(zeros(nPixels,nBatch));
avgMotion = single(zeros(nPixels,nBatch));
w = zeros(nBatch,1,'single');
wd = w;
parStartTime = clock;
if useParallel
    parfor i=1:nBatch
        fprintf('Starting batch #%d of %d...\n',i,nBatch)
        [avgFrame(:,i),avgMotion(:,i),w(i),wd(i)] = ...
            avgROI(vobj,indices(i),roi,nt,ns);
        fprintf('Rreading batch #%d of %d completed\n',i,nBatch)
    end
else
    for i=1:nBatch
        fprintf('Starting batch #%d of %d...\n',i,nBatch)
        [avgFrame(:,i),avgMotion(:,i),w(i),wd(i)] = ...
            avgROI(vobj,indices(i),roi,nt,ns);
        fprintf('Rreading batch #%d of %d completed\n',i,nBatch)
    end
end
parEndTime = clock;
fprintf('Averaging video frames from %s completed.\n',...
    videoFile{1})
fprintf('It took %6.2f.\n',etime(parEndTime,parStartTime))

% weighted mean by batch size
vSVD.avgFrame = (avgFrame*w)/sum(w);
vSVD.avgMotion = (avgMotion*wd)/sum(wd);

clear avgFrame avgMotion

% SVD on video ROI
uMot = cell(1,nBatch);
parStartTime = clock;
if useParallel
    parfor i=1:nBatch
        fprintf('Getting SV for batch #%d of %d...\n',i,nBatch)
        uMot{i} = doSVD(vobj,indices(i),roi,vSVD.avgMotion,nt,ns);
        fprintf('SV for batch #%d of %d completed\n',i,nBatch)
    end
else
    for i=1:nBatch
        fprintf('Getting SV for batch #%d of %d...\n',i,nBatch)
        uMot{i} = doSVD(vobj,indices(i),roi,vSVD.avgMotion,nt,ns);
        fprintf('SV for batch #%d of %d completed\n',i,nBatch)
    end
end
parEndTime = clock;
fprintf('SVD for session %s  completed.\n',videoFile{1})
fprintf('It took %6.2f.\n',etime(parEndTime,parStartTime))

ncomps = 1000;
idx = [indices.nFrames] == batchSize;
uMot = cell2mat(uMot(idx));

tic
[u,s,v] = svdecon(uMot);
toc
uMotMask = u(:,1:min(ncomps,size(u,2)));
vSVD.uMEraw = uMotMask;
vSVD.uME = normc(uMotMask);
SVD.uME = u;
SVD.sME = diag(s);
SVD.vME = v;
SVD.uMot = uMot;

clear u s v uMot uMotMask

% Projection of SVs on video frames to calculate tempral profiles of SV
% weights
motSVD = cell(nBatch,1);
me = cell(1,nBatch);
parStartTime = clock;
if useParallel
    parfor i=1:nBatch
        fprintf('Projecting on SV for batch #%d of %d...\n',i,nBatch)
        [motSVD{i}, me{i}] =  ...
            projectSVD(vobj,indices(i),roi,...
            vSVD.avgMotion,vSVD.uME,nt,ns);
        fprintf('Projection for batch #%d of %d completed\n',i,nBatch)
    end
else
    for i=1:nBatch
        fprintf('Projecting on SV for batch #%d of %d...\n',i,nBatch)
        [motSVD{i}, me{i}]  =  ...
            projectSVD(vobj,indices(i),roi,...
            vSVD.avgMotion,vSVD.uME,nt,ns);
        fprintf('Projection for batch #%d of %d completed\n',i,nBatch)
    end
end
parEndTime = clock;
fprintf('Projection session %s  completed.\n',videoFile{1})
fprintf('It took %6.2f.\n',etime(parEndTime,parStartTime))

nFramesDown = zeros(1,nVideoFiles);
if nVideoFiles == 1
    ME = cell2mat(me);
    vSVD.svdME = cell2mat(motSVD);
    nFramesDown = sum(w);
else
    vSVD.svdME = cell(nVideoFiles,1);
    ME = cell(1,nVideoFiles);
    for i=1:nVideoFiles
        idx = [indices.i] == i;
        ME{i} = cell2mat(me(idx));
        vSVD.svdME{i} = cell2mat(motSVD(idx));
        nFramesDown(i) = sum(w(idx));
    end
end
vSVD.vidHeightWidth = [vobj{1}.Height,vobj{1}.Width];
vSVD.ROI = roi;
vSVD.frameRate = vobj{1}.FrameRate;
vSVD.downsampleT = nt;
vSVD.downsampleS = ns;
vSVD.numFrames = nFrames;
vSVD.numDsamplFrames = nFramesDown;
vSVD.fileNm = videoFile;

end


function [avgFrame,avgMotion,w,wd] = avgROI(vobj,indices,mask,nt,ns)
roi = readOneBatch(vobj,indices,mask,nt,ns);
d = abs(diff(roi,1,2));

avgFrame = mean(roi,2);
avgMotion = mean(d,2);

w = size(roi,2);
wd = size(d,2);
end


function roi = readOneBatch(vobj,indices,mask,nt,ns)
nh = floor(vobj{1}.Height/ns);
nw = floor(vobj{1}.Width/ns);
idxh = 1:floor(vobj{1}.Height/ns)*ns;
idxw = 1:floor(vobj{1}.Width/ns)*ns;

if rem(indices.nFrames,nt) ~= 0
    error('Batch size should be a multiple of %d',nt)
end
ncol = indices.nFrames/nt;

nPixels = sum(mask(:));

roi = zeros(nPixels,ncol,'single');
v = read(vobj{indices.i},indices.index);
v = single(squeeze(v(:,:,1,:)));

if ns > 1
    for i=1:size(v,3)
        a = v(:,:,i);
        a = squeeze(mean(mean(reshape(a(idxh,idxw,1),...
            ns, nh, ns, nw), 1),3));
        roi(:,i) = a(mask);
    end
else
    for i=1:size(v,3)
        a = v(:,:,i);
        roi(:,i) = a(mask);
    end
end

if nt > 1
    roi = squeeze(mean(reshape(roi,size(roi,1),nt,ncol),2));
end
end


function uMot = doSVD(vobj,indices,mask,avgMotion,nt,ns)
ncomps = 200;

roi = readOneBatch(vobj,indices,mask,nt,ns);

imot = bsxfun(@minus,abs(diff(roi,1,2)),avgMotion);
imot = gpuArray(imot);

u = svdecon(imot);
u = gather_try(u);

uMot = u(:,1:min([ncomps,size(u,2)]));
end


function [motSVD,me] = projectSVD(vobj,indices,mask,avgMotion,uMotMask,nt,ns)
if indices.index(2) <= vobj{indices.i}.NumFrames - nt
    indices.index(2) = indices.index(2) + nt;
    indices.nFrames = indices.nFrames + nt;
end

roi = readOneBatch(vobj,indices,mask,nt,ns);
d = abs(diff(roi,1,2));
motSVD = (d-avgMotion)'*uMotMask;
me = mean(d,1);
end



