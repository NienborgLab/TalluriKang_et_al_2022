% for face video
videoFile = {'faceVideo.avi'};
roiFile = {'faceROI.mat'};
[face.vSVD,face.SVD,face.ME] = motionSVD(videoFile, roiFile);

% for body video - illustrates downsampling
videoFile = {'bodyVideo.avi'};
roiFile = {'bodyROI.mat'};
[body.vSVD,body.SVD,body.ME] = motionSVD(videoFile, roiFile,3,2);
