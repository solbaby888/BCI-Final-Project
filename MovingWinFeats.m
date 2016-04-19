function feat = MovingWinFeats(x, fs, winLen, winDisp, featFn)

%====READ ME==============================================================
%{
    Calculates  feature of a sample with a given sliding window
    winLen and winDisp must be seconds fs must be in Hz 
    xLen must be in sample
%}
%=========================================================================


% Calculate featfn Feature

% input values
xLen      = length(x);
NoW       = 1+floor((xLen-fs*winLen)/(winDisp*fs));
win       = fs*winLen; % window (samples)
win_d     = fs*winDisp;

% initial values
first     = 1;
last      = win;

% intializing
feat = [];

% calculate LL for each sliding window and store in array
for i = 1:NoW
    SW         = x(first: last); % sliding window
    feat(i,:)    = featFn(SW, fs);     % calculate feature
    first      = first + win_d;
    last       = last  + win_d;
end;
    
    
