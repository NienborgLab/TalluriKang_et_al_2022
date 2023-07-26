function [motion_idx] = findMotion(vel,fps,cri,tol)

% simplified version of Incheol's function 'findMotion.m'

%input:
% vel       vector of motion energy or velocity
% fps       sampling frequency in Hz
% cri       threshold above which we define motion to have occurred
% tol       tolerance (in sec) to combine motions occuring within this
%           tolerance interval
%
% output
% motion_idx  n-by-2 matrix of on (col 1) and off (col 2) indices for
%            motion epochs

if nargin == 3
    tol = 1;
end

d = diff(vel >= cri);

idx_on = find(d==1)+1;
idx_off = find(d==-1);

% adjust for unequal number of indices
if length(idx_on)>length(idx_off)
    if idx_on(2)<idx_off(1)
        idx_on = idx_on(2:end);
    else
        idx_on = idx_on(1:end-1);
    end
end
if length(idx_off)>length(idx_on)
    if idx_off(1)<idx_on(1)
        idx_off = idx_off(2:end);
    else
        idx_off = idx_off(1:end-1);
    end
end
motion_idx = [idx_on(:),idx_off(:)];

% inter-motion intervals
% combine motions occurring within tolerance interval
d = (motion_idx(2:end,1) - motion_idx(1:end-1,2)) / fps;
idx = d < tol;

motion_idx([false;idx],1) = NaN;
motion_idx([idx;false],2) = NaN;
idx = all(isnan(motion_idx),2);
motion_idx(idx,:) = [];

idx = find(isnan(motion_idx(:,1)));
motion_idx(idx-1,2) = motion_idx(idx,2);
motion_idx(idx,:) = [];


