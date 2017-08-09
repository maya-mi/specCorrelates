
load('Sun/timesSun')


load('ironA.mat')

load Sun/processedSun.mat
load Sun/fitResultsSun
r = mean(reduced);
redFilter = find(r < 2);

goodLines = ironA(redFilter);
dlmwrite('Sun/goodLines.csv', goodLines', 'precision', 7)

for k = 1:length(goodLines)
    fig = figure;
    set(fig, 'Position', [400 0 600 600])
    i = redFilter(k);
    subplot(2, 2, 1)
    scatter(mean(squeeze(wavelengths(:, i, :))), mean(squeeze(normOrders(:, i, :))))
    title(ironA(i))
    subplot(2, 2, 2)
    scatter(uniqueNights, squeeze(f(:, i, 3)))
    title('Center Offset')
    subplot(2, 2, 3)
    scatter(uniqueNights, squeeze(f(:, i, 4)))
    title('Width (std)')
    subplot(2, 2, 4)
    scatter(uniqueNights, squeeze(f(:, i, 2)./ f(:, i, 1)))
    title('Relative Depth')
end
    
error('Update Sun/goodLines.csv 2nd column with 1/0 good/bad filter')





