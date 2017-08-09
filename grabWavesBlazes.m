load('ironLineCoords.mat')

numlines = length(coords);

import matlab.io.*
files = dir('**/*wave_A.fits');
numfiles = length(files);

widthRange = 15;
allWaves = zeros(numfiles, numlines, 2*widthRange + 1);
waveTimes = zeros(numfiles, 1);

for k = 1: numfiles
  fName = strcat(files(k).folder, '/', files(k).name);
  f = fits.openFile(fName);
  try
      [t, ~] = fits.readKey(f, 'DATE-OBS');
      waveTimes(k) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
      fptr = fitsread(fName);
      for j = 1: numlines
          allWaves(k, j, :) = fptr(coords(j, 1), coords(j, 2) - widthRange: coords(j, 2) + widthRange);
      end
  catch
      waveTimes(k) = nan;
  end
  
  fits.closeFile(f)
end
catchFilter = ~isnan(waveTimes);
waveTimes = waveTimes(catchFilter);
allWaves = allWaves(catchFilter, :, :);

save('Sun/ironWaves.mat', 'allWaves', 'waveTimes', '-v7.3')



files = dir('**/*blaze_A.fits');
numfiles = length(files);
blazeTimes = zeros(numfiles, 1);
allBlazes = zeros(numfiles, numlines, 2*widthRange + 1);

for k = 1: numfiles
  fName = strcat(files(k).folder, '/', files(k).name);
  f = fits.openFile(fName);
  try
      [t, ~] = fits.readKey(f, 'DATE-OBS');
      blazeTimes(k) = datenum(datetime(t, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
      fptr = fitsread(fName);
      for j = 1: numlines
          allBlazes(k, j, :) = fptr(coords(j, 1), coords(j, 2) - widthRange: coords(j, 2) + widthRange);
      end
  catch
      blazeTimes(k) = nan;
  end
  fits.closeFile(f)
end

catchFilter = ~isnan(blazeTimes);
blazeTimes = blazeTimes(catchFilter);
allBlazes = allBlazes(catchFilter, :, :);

save('Sun/ironBlazes.mat', 'allBlazes', 'blazeTimes', '-v7.3')


