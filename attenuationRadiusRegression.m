function estimatedAttenutation = attenuationRadiusRegression(manualSelectFlag,polyFitN)
%% attenuationRadiusRegression.m
% this function fits a polynomial to the function which maps pullback distance to attenutation
% (signal x).  in the future, once particular geometrical data is available, it can be generalized
% to regress some function (not necessarily polynomial) between radial distance from "plaque" to
% source/detector and attenutation

%% load data
loadDistancePhantom

%% manual selection of training set
switch manualSelectFlag
    case true
        TUBEANGLE = 7;  %degrees, replace hard code later
        TUBEMINDISTANCE = .5;   %mm, from thesis
        plot(signal);
        [idx,attenuationPlaque] = getpts;
        idx = floor(idx);
        pullBackDistance = [0:ltr/(length(signal)-1):ltr]';
        radius = TUBEMINDISTANCE + tand(TUBEANGLE)*pullBackDistance;
        radiusPlaque = radius(idx);
        %save('attenuationAndRadiiOfPlaqueEvents','attenuationPlaque','radiusPlaque','idx','pullBackDistance','radius');
    case false
        load('attenuationAndRadiiOfPlaqueEvents.mat');
end

%% polynomial fitting
polyAttenuationRadius = polyfit(1./radiusPlaque,attenuationPlaque,polyFitN);
estimatedAttenutation = polyval(polyAttenuationRadius,1./radius);

%% check poly fit (graph)
figure;
plot(radius,estimatedAttenutation,'r');
xlabel('Radius (mm)');
ylabel('Attenuation');
hold on;
plot(radius,signal,'b');
scatter(radiusPlaque,attenuationPlaque,'b*')
legend({'estimated Attenutation w/ plaque', 'measured attenuation (possibly no plaque)'});

end
