load('ironLineCoords.mat')
numlines = length(coords);
c = 299792.458; %speed of light in km/sec
import matlab.io.*
files = dir('**/*e2ds_A.fits');
numfiles = length(files);
numObs = length(directory);
obsTimes = zeros(numObs, 1);
berv = zeros(numObs, 1); 
airmass = zeros(numObs, 1);
julian = zeros(numObs, 1);
exptime = zeros(numObs, 1);
widthRange = 15;
pix = zeros(numObs, numlines, 2*widthRange + 1);
allLines = zeros(numObs, numlines, 2*widthRange + 1);
allWaves = zeros(numObs, numlines, 2*widthRange + 1);
allBlazes = zeros(numObs, numlines, 2*widthRange + 1);
waveTimes = zeros(numObs, 1);
blazeTimes = zeros(numObs, 1);
currentFolder = 'notAssigned';

for i = 1: numObs
    %Use pre-made directory of obs to only open relevant files
    k = directory(i);
    fName = strcat(files(k).folder, '/', files(k).name);
    f = fits.openFile(fName);
    val = fits.readKey(f, 'OBS-TYPE');
    if strcmp(val,'SCIENCE') && strcmp(fits.readKey(f, 'HIERARCH TNG OBS TARG NAME'), starName)
      t = fits.readKey(f, 'DATE-OBS');
      obsTimes(i) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
      berv(i) = fits.readKeyDbl(f, 'HIERARCH TNG DRS BERV');
      airmass(i) = fits.readKeyDbl(f, 'AIRMASS');
      julian(i) = fits.readKeyDbl(f, 'MJD-OBS');
      exptime(i) = fits.readKeyDbl(f, 'EXPTIME');
      exp = fitsread(fName);
      %grab wave, blaze files for new folder, if changed
      if ~strcmp(files(k).folder, currentFolder)
          currentFolder = files(k).folder;
          dirWv = dir(strcat(files(k).folder, '/*wave_A.fits'));
          dirBlz = dir(strcat(files(k).folder, '/*blaze_A.fits'));
          waveTimesObs = zeros(length(dirWv), 1);
          blazeTimesObs = zeros(length(dirBlz), 1);
          for w = 1:length(dirWv)
              wName = strcat(files(k).folder, '/', dirWv(w).name);
              wave = fits.openFile(wName);
              t = fits.readKey(wave, 'DATE-OBS');
              waveTimesObs(w) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
              fits.closeFile(wave)
          end
          for b = 1:length(dirBlz)
              bName = strcat(files(k).folder, '/', dirBlz(b).name);
              blaze = fits.openFile(bName);
              t = fits.readKey(blaze, 'DATE-OBS');
              blazeTimesObs(b) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
              fits.closeFile(blaze)
          end
      end
      [~, W] = min(abs(waveTimesObs - obsTimes(i)));
      [~, B] = min(abs(blazeTimesObs - obsTimes(i)));
      obsWave = fitsread(strcat(files(k).folder, '/', dirWv(W).name));
      obsBlaze = fitsread(strcat(files(k).folder, '/', dirBlz(B).name));
      waveTimes(i) = waveTimesObs(W);
      blazeTimes(i) = blazeTimesObs(B);
      for j = 1: numlines
          try
          bervEff = berv(i) - starOffset;
          lambdaNew = ironA(j)*(1 - bervEff/c)/ sqrt(1 - (bervEff/c)^2);
          waveStep = mean(diff(obsWave(coords(j, 1), coords(j, 2) - widthRange:coords(j, 2)+widthRange)));
          deltaPix = round((lambdaNew - ironA(j)) / waveStep);
          rAdj = coords(j, 2) + deltaPix;
          grabPix = rAdj - widthRange: rAdj + widthRange;
          pix(i, j, :) = grabPix;
          allLines(i, j, :) = exp(coords(j, 1), grabPix);
          allWaves(i, j, :) = obsWave(coords(j, 1), grabPix);
          allBlazes(i, j, :) = obsBlaze(coords(j, 1), grabPix);
          catch
          end
      end
    end
  fits.closeFile(f)  
end


save (strcat(starName, '/WaveBlazes', starName, '.mat'),'allWaves', 'allBlazes', 'waveTimes', 'blazeTimes', 'pix', 'berv','-v7.3')
save (strcat(starName, '/ironLines', starName, '.mat'), 'allLines', 'julian', '-v7.3')
save (strcat(starName, '/headerExtras', starName, '.mat'), 'airmass', 'exptime', 'obsTimes')
   