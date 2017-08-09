nOrders = 69;
load('ironA.mat') %file containing vector of iron line wavelengths in Angstroms
load('avgWaves.mat')

%making order bins, discretizing into bins
bins = [min(avgWaves, [], 2); max(avgWaves(end, :))];
d = discretize(ironA, bins);
% hasLines = ismember(1:nOrders, d);
coords = [];
for i = 1:nOrders
    orderLines = ironA(d == i);
    for j = 1:length(orderLines)
        [~, col] = min(abs(orderLines(j) - avgWaves(i, :)));
        coords = [coords; i , col];
    end
end
ironA = ironA(coords(:, 2) > 30 & coords(:,2) < 4066, :);
coords = coords(coords(:, 2) > 30 & coords(:,2) < 4066, :);

save('ironLineCoords.mat', 'coords', 'ironA');
  