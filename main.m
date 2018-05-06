%% distancePhantomAnalysis
% this script implements a matched filter detection of plaque events in the distance phantom.
% specifically, it detects a deterministic gaussian (supported by Mallas thesis on the assumption of
% a gaussian point spread function from plaque source across angle) whose width is fit to some
% manually selected training set of plaque events in signal (see findOptimalWidth) and whose
% amplitude varies with radius (here it is implemented as varying with pullback distance as a proxy,
% though given particular geometry it can be extended later).  We fit a polynomial to the pullback
% distance to attenuation points from plaque events manually selected by a user (see
% attenuationRadiusRegression)

%This implementation is not computationally efficient!

%% parameters
clearvars; clc;
PFA = 10^-1;

%% compute attenuation dependance on pullback distance
%turns manual select on or off (off loads previously selected points)
manualSelectFlag = false;

%degree of polynomial to fit to distance vs attenutation of plaque
polyFitN = 2;

estimatedAttenutation = attenuationRadiusRegression(manualSelectFlag,polyFitN);

%% find best fit width of gaussian to data
% examine only the firstN known plaque events, as these are the most certain
onlyFirstN = 10;

% showWhich is a vector of indexes from 1 to onlyFirstN which determines which signals (and their
% fits) are graphed at the end of this function.  leaving it empty turns graphing off
showWhich = [1:4];

% templateSignalWidth is the length (in mm) of signal which we convolve in using the matched filter
templateSignalStats.widthMM = .1;

[templateSignalStats] = findOptimalWidth(onlyFirstN,showWhich,templateSignalStats);

%% load data
loadDistancePhantom;
load('attenuationAndRadiiOfPlaqueEvents.mat');

%% compute noiseSigma based on last 1/12-th of data set
startIdx = round(length(signal)*9/12);
noiseSigma = std(signal(startIdx:length(signal)));

%% compute test statistic, T, threshold and probability (assuming uniform prior) that signal is from plaque event
signalPower = @(x)(sum(x.^2));

%functions taken from appendix 2c (page 51) of Steven Kay Detection Theory (Vol 2)
Q = @(x)(1-normcdf(x));
Qinv = @(x)(-norminv(x,0,1));

dummyDistance = pullBackDistance(1:templateSignalStats.widthSamples);
dummyDistance = dummyDistance - mean(dummyDistance); %just to give scale, estimated mu wil be wrong

%compute template and its power
template = normpdf(dummyDistance,0,templateSignalStats.sigma);
template = template/max(template);
templatePower = signalPower(template);

for pullBackIdx = templateSignalStats.halfWidth+1:length(signal)-templateSignalStats.halfWidth

    templatePowerGivenRadius(pullBackIdx) = templatePower*estimatedAttenutation(pullBackIdx)^2;
    
    %compute test statistic
    scale(pullBackIdx) = sqrt(templatePowerGivenRadius(pullBackIdx)*noiseSigma^2);
    startIdx = pullBackIdx-templateSignalStats.halfWidth;
    endIdx = pullBackIdx+templateSignalStats.halfWidth;
    T(pullBackIdx) = signal(startIdx:endIdx)*template/scale(pullBackIdx);    %no need to flip, gauss is symmetric
    
    %compute threshold
    gamma(pullBackIdx) = scale(pullBackIdx)*Q(PFA);
    
    %compute Pd
    Pd(pullBackIdx) = Q(Qinv(PFA)-sqrt(signalPower(template)/noiseSigma^2));
    
    %compute prob(plaque|attenuation) assuming a uniform prior (incomplete)
    
end

Pd = Pd(Pd~=0);     %this is ugly but the indexing is restricted in parfor loops...

%% graphs

% detected vs not detected (original signal) plot red - plaque event, blue - no plaque event
plaquePresent = T > gamma;

plaqueSignal = signal;  plaqueSignal(~plaquePresent) = NaN;
noPlaqueSignal = signal;  noPlaqueSignal(plaquePresent) = NaN;
figure;
plot(radius,plaqueSignal,'r');
xlabel('radius (mm)');
ylabel('Attenuation');
hold on;
plot(radius,noPlaqueSignal,'b');
legend({'plaque detected', 'plaque not-detected'});

% detection as an image
figure
fullPlaquePresent = [nan(1,100),plaquePresent,nan(1,100)];
detectionImage = reshape(fullPlaquePresent,2000,60)';
detectionImage = detectionImage(6:size(detectionImage,1),:);    %first 5 rows are startup data ...
imagesc(detectionImage);
set(gca,'XTickLabel',[36:36:360]);
xlabel('Angle (deg)');
set(gca,'YTickLabel',[min(radius):floor((max(radius)-min(radius))/10*100)/100 :max(radius)]);
ylabel('Target Distance (mm)');


% PD and pullback distance (not working)
figure;
clippedPullBackDistance = pullBackDistance(templateSignalStats.halfWidth+1:length(signal)-templateSignalStats.halfWidth);
plot(clippedPullBackDistance,Pd);
xlabel('Pullback Distance (mm)');
ylabel('Probability of detection');