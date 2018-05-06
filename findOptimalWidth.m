function [templateSignalStats, dummyDistance] = findOptimalWidth(onlyFirstN,showWhich,templateSignalStats)
%% findOptimalWidth
% this script finds the optimal width for the gaussian template signal by choosing the width
% corresponding to the max point of convolution of template gaussians (of varying amplitude and
% widths) with a segment of the signal which is known to have a plaque response.  (signal is
% manually segmented first)

%% load data
loadDistancePhantom;
load('attenuationAndRadiiOfPlaqueEvents.mat');     %manually selected attenuation/distance pts from known plaque

%% derived parameters
templateSignalStats.samplesPerMM = length(signal)/ltr;              %ltr = linear travel length
templateSignalStats.widthSamples = round(templateSignalStats.samplesPerMM*templateSignalStats.widthMM);
templateSignalStats.widthSamples = templateSignalStats.widthSamples  + ~mod(templateSignalStats.widthSamples,2);   %make odd
templateSignalStats.halfWidth = floor(templateSignalStats.widthSamples/2);

%% build known plaque test cases, compute their optimal gaussian fit
dummyDistance = pullBackDistance(1:templateSignalStats.widthSamples);     %just to give scale, estimated mu wil be wrong
for testCaseIdx = onlyFirstN:-1:1
    % segment out known plaque test case
    startIdx = idx(testCaseIdx)-templateSignalStats.halfWidth;
    endIdx = idx(testCaseIdx)+templateSignalStats.halfWidth;
    testSignal(testCaseIdx,:) = signal(startIdx:endIdx);
    
    %compute optimal gaussian fit
    [sigma(testCaseIdx),mu(testCaseIdx),amplitude(testCaseIdx)] = mygaussfit(dummyDistance,testSignal(testCaseIdx,:));
end

%% choose width corresponding to max mean across segmented signals
templateSignalStats.sigma = mean(sigma);

%% plot as sanity check
figure;
for testCaseIdx = 1:length(showWhich)
    plot(dummyDistance,testSignal(testCaseIdx,:),'b');
    hold on;
    pdf = normpdf(dummyDistance,mu(testCaseIdx),sigma(testCaseIdx));
    pdf = pdf/max(pdf)*amplitude(testCaseIdx); % rescale
    plot(dummyDistance,pdf,'r');
    xlabel('Distance (mm), note all signals shifted');
    ylabel('Attenuation');
end

end