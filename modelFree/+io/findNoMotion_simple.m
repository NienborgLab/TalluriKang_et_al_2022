function [motion_idx] = findNoMotion(vel,fps,cri,tol)

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

idx_on = find(d==-1)+1; % period below motion threshold
idx_off = find(d==1);

if vel(1)<cri  % we are starting below the citerion
    idx_off = idx_off(2:end);
end

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


% inter-no-motion intervals
% remove no-motion that is within tolerance of no-motion onset
idx_on = idx_on + round(tol*fps);
idx_off = idx_off-round(tol*fps);

incl = find(idx_on <length(vel) & idx_off>1);
idx_on = idx_on(incl);
idx_off = idx_off(incl);

d = idx_off-idx_on;

idx = d>0;


if size(idx_on,2) >size(idx_on,1)
    idx_on = idx_on';
    idx_off = idx_off';
end

motion_idx = [idx_on(idx,1),idx_off(idx,1)];

% % combine motions occurring within tolerance interval
% d = (motion_idx(2:end,1) - motion_idx(1:end-1,2)) / fps;
% idx = d < tol;
% 
% motion_idx([false;idx],1) = NaN;
% motion_idx([idx;false],2) = NaN;
% idx = all(isnan(motion_idx),2);
% motion_idx(idx,:) = [];
% 
% idx = find(isnan(motion_idx(:,1)));
% motion_idx(idx-1,2) = motion_idx(idx,2);
% motion_idx(idx,:) = [];
% 

