load('ironLineCoords.mat')
numlines = length(coords);

load targets
directory = find(strcmp(target, 'Sun'));
numObs = length(directory);

import matlab.io.*
files = dir('**/*e2ds_A.fits');

obsTimes = zeros(numObs, 1);
berv = zeros(numObs, 1); 
airmass = zeros(numObs, 1);
julian = zeros(numObs, 1);
exptime = zeros(numObs, 1);
widthRange = 15;
allLines = zeros(numObs, numlines, 2*widthRange + 1);


for i = 1: numObs
    k = directory(i);
    fName = strcat(files(k).folder, '/', files(k).name);
    f = fits.openFile(fName);
    val = fits.readKey(f, 'OBS-TYPE');
    tar = fits.readKey(f, 'HIERARCH TNG OBS TARG NAME');
    if strcmp(val,'SCIENCE') && strcmp(tar, 'Sun')
      t = fits.readKey(f, 'DATE-OBS');
      obsTimes(i) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
      berv(i) = fits.readKeyDbl(f, 'HIERARCH TNG DRS BERV');
      airmass(i) = fits.readKeyDbl(f, 'AIRMASS');
      julian(i) = fits.readKeyDbl(f, 'MJD-OBS');
      exptime(i) = fits.readKeyDbl(f, 'EXPTIME');
      fptr = fitsread(fName);
      for j = 1: numlines
          allLines(i, j, :) = fptr(coords(j, 1), coords(j, 2) - widthRange: coords(j, 2) + widthRange);
      end
    end
  fits.closeFile(f)  
end


save('Sun/ironLinesSun.mat', 'allLines', 'julian', 'exptime', '-v7.3')
save('Sun/headerExtrasSun.mat', 'airmass', 'obsTimes', 'berv', '-v7.3')
