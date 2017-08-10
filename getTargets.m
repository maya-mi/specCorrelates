load('ironLineCoords.mat')
numlines = length(coords);
c = 299792.458; %speed of light in km/sec
import matlab.io.*
files = dir('**/*e2ds_A.fits');
numfiles = length(files);
target = cell(numfiles, 1);
julian = zeros(numfiles, 1);

display(strcat(num2str(numfiles), ' total files'))

for k = 1: numfiles
  if mod(k, 1000) == 0
      display(strcat(num2str(k/numfiles*100), '%'))
  end
  try
    fName = strcat(files(k).folder, '/', files(k).name);
    f = fits.openFile(fName);
    julian(k) = fits.readKeyDbl(f, 'MJD-OBS');
    val = fits.readKey(f, 'OBS-TYPE');
    if strcmp(val,'SCIENCE')
        target{k} = fits.readKey(f, 'HIERARCH TNG OBS TARG NAME');
    else
        target{k} = 'notScience';
    end  
  catch
      target{k} = 'BROKEN';
  end
  try
    fits.closeFile(f) 
  catch
  end
end

save targets.mat target julian
   