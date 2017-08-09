nOrders = 69;
cs = 4096;
load(strcat(starName, '/WaveBlazes', starName, '.mat'))
load(strcat(starName, '/IronLines', starName, '.mat'))
load ironA.mat

[nObs, nLines, nPix] = size(allLines);

%Prepping for the dark magic of blaze scaling
slope = -0.32557502060181/4000;
intercept = 1 - slope * 4095/2;
blazeAdjScalars = slope*pix +intercept;
blazeAdj = allBlazes .* blazeAdjScalars;
sortedScales = sort(allLines ./ blazeAdj, 3);
%Pick ~top 3rd of ratios-- as a line is about 15 points wide, this serves
%as a good proxy for local continuum value if used with numPix = 31
blazeScales = mean(sortedScales(:, :, round(nPix * 2 / 3):nPix), 3);
blazeFinal = blazeScales .* blazeAdj;
normOrders = allLines ./ blazeFinal;

% Propagating error:
blazeError = 1.0000e-03 * blazeFinal; 
propError = sqrt(1 ./ allLines + (blazeError ./ blazeFinal) .^2) .* normOrders;
save(strcat(starName, '/propError', starName, '.mat'), 'propError', 'blazeError', '-v7.3')


%JPL correction
c = 299792.458; %speed of light in km/sec
bervEff = berv - starOffset;
wavelengths = allWaves .* (1 + bervEff / c) ./ sqrt(1 - (bervEff/c).^2);

save(strcat(starName, '/processed', starName, '.mat'), 'normOrders', 'ironA', 'wavelengths', '-v7.3')

obsTimes = datetime(julian + 2.4e+6, 'ConvertFrom', 'julianDate');
obsNights = dateshift(obsTimes - hours(12), 'start', 'day');
uniqueNights = unique(obsNights);
save(strcat(starName, '/times', starName, '.mat'), 'obsTimes', 'obsNights', 'uniqueNights');
