nOrders = 69;
cs = 4096;

load 'Sun/ironLinesSun.mat'
obsTimes = datetime(julian + 2.4e+6, 'ConvertFrom', 'julianDate');

%Cloud cuts: 
T = readtable('Sun/RVs3584099468');
cloudOK = datetime(table2array(T(:, 1)), 'ConvertFrom', 'juliandate');

%Finds appt times w/ tolerance of 3 sec
load Sun/headerExtrasSun
%midExp = obsTimes + seconds(exptime/2);
midExp = obsTimes + minutes(2.5); %enforcing 5-min exposures
noClouds = datefind(cloudOK, midExp, seconds(3)); 
obsTimes = datetime(obsTimes(noClouds), 'ConvertFrom', 'datenum');

%The names here are funny for the sun, but consistent with star routine
obsNights = dateshift(obsTimes, 'start', 'day');
uniqueNights = unique(obsNights);

%for SDO comparisons:
load('Sun/SDOFIt_Parma_2015_2017')
load('Sun/SDOFIt_Parma_2015_2017_GeomFactor.mat');
load('Sun/Bflux2015_2017_Data_RunAt_31-Jul-2017_22_01_29');
sdoTimes = datetime(DaTSDO, 'ConvertFrom', 'datenum');
sdoNights = dateshift(sdoTimes, 'start', 'day');
[hasSDO, ~] = ismember(uniqueNights, sdoNights);

uniqueNights = uniqueNights(hasSDO);
sdoFilter = ismember(obsNights, uniqueNights);
obsNights = obsNights(sdoFilter);
obsTimes = obsTimes(sdoFilter);

%Mandating 5 exposures per day: 
edges = [uniqueNights; max(uniqueNights) + days(1)];
[N, ~] = histcounts(obsNights, edges);
uniqueNights = uniqueNights(N > 5);
numsFilter = ismember(obsNights, uniqueNights);
obsNights = obsNights(numsFilter);
obsTimes = obsTimes(numsFilter);

save('Sun/timesSun.mat', 'obsTimes', 'obsNights', 'uniqueNights');


%Applying same filters to allLines to keep times lined up:
allLines = allLines(noClouds, :, :);
allLines = allLines(sdoFilter, :, :);
allLines = allLines(numsFilter, :, :);

%Finding SDO endpoints for comparison: Note that days start at midnight
nDays = length(uniqueNights);
filling = zeros(nDays, 1);
plague = zeros(nDays, 1); 
spot = zeros(nDays, 1); 
bFlux = zeros(nDays, 1); 
vs = zeros(nDays, 3);
vsLabels = {'vCon', 'vPhoto', 'vQuiet'};
for k = 1: nDays
    sameDay = sdoNights == uniqueNights(k);
    filling(k) = mean(ActFrac(sameDay));
    plague(k) = mean(PlagFrac(sameDay));
    spot(k) = mean(SpotFrac(sameDay));
    bFlux(k) = mean(Bflux(sameDay));
    vs(k, :) = [mean(v_con(sameDay)), mean(v_phot(sameDay)), mean(v_quiet(sameDay))];
end

save('Sun/sdoComp.mat', 'filling', 'plague','spot', 'vs', 'vsLabels', 'bFlux');


%For RV correction
load('Sun/JPL_Horizons_1min2.mat')
jplTimes = datetime(DateTimeJPL, 'ConvertFrom', 'datenum');

%Scaling to time wave, blaze and JPL RV:
[nObs, nLines, nPix] = size(allLines);
blazeI = zeros(nObs, 1); %tracks indices for blazes
waveI = zeros(nObs, 1); 
rvI = zeros(nObs, 1);

%BLAZES AND WAVES STORED AS 3D ARRAYS blazes AND waves
load('Sun/ironBlazes.mat')
load('Sun/ironWaves.mat')
waveTimes = datetime(waveTimes, 'ConvertFrom', 'datenum');
blazeTimes = datetime(blazeTimes, 'ConvertFrom', 'datenum');
%Deal with any out-of-order times:
[waveTimes, wI] = sort(waveTimes);
[blazeTimes, bI] = sort(blazeTimes);
allWaves = allWaves(wI, :, :);
allBlazes = allBlazes(bI, :, :);

%Find location of nearest wave/blaze/rv time to each exposure
for k = 1:nObs
    [~, rvI(k)] = min(abs(obsTimes(k) - jplTimes));
    [~, blazeI(k)] = min(abs(obsTimes(k) - blazeTimes));
    [~, waveI(k)] = min(abs(obsTimes(k) - waveTimes));
end

%Find corresponding vals
blazes = allBlazes(blazeI, :, :);
waves = allWaves(waveI, :, :);
rvs = RV_JPL(rvI);


%Prepping for the dark magic of blaze scaling--

%Get col coord of points for linear blaze adjustment
load('ironLineCoords.mat')
pixNums = zeros(1, nLines, nPix);
for i = 1:nLines
    pixNums(1, i, :) = coords(i, 2) - 15:coords(i, 2) + 15;
end

slope = -0.32557502060181/4000;
intercept = 1 - slope * 4095/2;
blazeAdjScalars = (slope*pixNums +intercept);
blazeAdj = blazes .* blazeAdjScalars;
sortedScales = sort(allLines ./ blazeAdj, 3);
%Pick ~top 3rd of ratios-- as a line is about 15 points wide, this serves
%as a good proxy for local continuum value if used with numPix = 31
blazeScales = mean(sortedScales(:, :, round(nPix * 2 / 3):nPix), 3);
blazeFinal = blazeScales .* blazeAdj;
normOrders = allLines ./ blazeFinal;

% Propagating error:
blazeError = 1.0000e-03 * blazeFinal; 
propError = sqrt(1 ./ allLines + (blazeError ./ blazeFinal) .^2) .* normOrders;
save('Sun/propErrorSun.mat', 'propError', 'blazeError', '-v7.3')


%JPL correction
c = 299792.458; %speed of light in km/sec
wavelengths = waves .* (1 - rvs / c);

save('Sun/processedSun.mat', 'normOrders', 'ironA', 'wavelengths', '-v7.3')

%Give primary guess for sun fit
widths = .05 * ones(1, nLines);
save widths.mat widths
