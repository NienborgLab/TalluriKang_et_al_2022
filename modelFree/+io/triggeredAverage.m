function [xT,xTA,t] = triggeredAverage(trig,signal,tWin,fS,offset)
% function [xT,xTA,t] = triggeredAverage(trig,signal,tWin,fS,offset)
% 
% input:
% trig      vector of indices around which we want to align the signal
% signal    signal for which to compute a triggered average
% tWin      temporal window (in #of samples) for which to compute the xTA: -tWin(1):tWin(2)
% fS        sampling frequency (Hz) of the signal (default is 60)
% offset    offset in time (sec) to align the signal to the triggers
%
% output:
% xT        n-by-m matrix where n is the number of triggers included
% xTA       x triggered average signal 
% t         time vector for the triggered average, aligned on
%           trigger+offset

% history
% 10/21/21  hn: wrote it

% parse inputs
if nargin ==3
    fS          =60;
    offset      = 0;
elseif nargin ==4
    offset      =0;
end

% exclude triggers for which the temporal window exceeds the available
% signal
idx             = find(trig>tWin(1) & trig<= length(signal)-tWin(2));
trig            = trig(idx);

% compute x-triggered average
xT              = zeros(length(trig),sum(tWin)+1);
for n = 1:length(trig)
    xT(n,:)     = signal(trig(n)-tWin(1):trig(n)+tWin(2));
end
xTA             = mean(xT,1);

t               = (-tWin(1):tWin(2))./fS + offset;

