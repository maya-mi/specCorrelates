
rvShift = [0, -0.354, -59.381, -7.752, -18.83, 26.78, -98.35];
stars = {'Sun', 'HD127334', 'HD144579', 'HD62613', 'HD219134', 'HD185144', 'HD103095'};
save starInfo.mat rvShift stars
load targets
for starCounter = 1:length(stars)
    starName = stars{starCounter};
    directory = find(strcmp(target, starName));
    starOffset = rvShift(starCounter);
    mkdir(starName);
    if strcmp(starName, 'Sun')
        grabIronLinesSun
        grabWavesBlazes
        normalizeSun
    else
        grabHDLinesBerv
        normalizeHD
    end
    clearvars -except starCounter starName stars rvShift target
    firstFit
    if strcmp(starName, 'Sun')
        load Sun/widthsOffsetsSun
        save widths.mat widths
    end
    fitHD
    clearvars -except starCounter stars rvShift target
end