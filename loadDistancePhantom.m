%% loadDistancePhantom.m

fileName = 'distance1.bin';
file = which(fileName);
[fid,~]=fopen(file,'rb','b');

params= fread(fid, 12,'double');
fs = params(1); % sampling frequency
ppr = params(3); % pixels per round
rpm = params(4); % rounds per minute
fdur = params(2);  % frame duration
step = params(6); % step size
numfr = params(7); % number of frames
ltr = params(5);  % linear travel range

% number of samples in one round
spr = 60*fs/rpm; 

signal = zeros(1, fs*fdur*numfr);
for t=1:numfr
    % 2*fs because we have data for the fluorescent and the intrinsic
    [A,~] = fread(fid, fs*2*fdur,'double');
    
    % fluorescent signal 
    % minus is for inverting the values to improve the image look
    A1 = -A(1:fs*fdur); 
    
    % stores cumulatively each frame's data in a row vector
    signal(fs*fdur*(t-1)+1 : fs*fdur*t) = A1';
end

%% zero mean the signal
signal = signal - mean(signal);