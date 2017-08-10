load Sun/sdoComp
load Sun/timesSun
load ironA
load Sun/lineFilter
load Sun/fitResultsSun.mat

[dCor, dSeries, dI, dSigns] = optimizeComp(squeeze(f(:, lineFilter, 2)./f(:, lineFilter, 1)), vs(:, 1));
[cCor, cSeries, cI, cSigns] = optimizeComp(squeeze(f(:, lineFilter, 3)), vs(:, 1));
[wCor, wSeries, wI, wSigns] = optimizeComp(abs(squeeze(f(:, lineFilter, 4))), vs(:, 1));

[d, dsc, ~] = normSeries(dSeries);
[c, csc, ~] = normSeries(cSeries);
[w, wsc, ~] = normSeries(wSeries);

[~, ~] = findWeightsCartesian(d, c, w, vs(:, 1));
title('Sun:vCon Weighting', 'Color', 'w', 'FontSize', 32)
[weights, bestCorr] = findWeightsParam(d, c, w, vs(:, 1));
title('Sun:vCon Weighting', 'Color', 'w', 'FontSize', 32)

T = readtable('Sun/RVs3584099468.txt');
sIDays = dateshift(datetime(T.obsJDcut, 'ConvertFrom', 'juliandate'), 'start', 'day');
nDays = length(uniqueNights);
sI = zeros(nDays, 1);
rvs = zeros(nDays, 1);
for D = 1:nDays
    sI(D) = mean(T.Sindexcut(sIDays == uniqueNights(D)));
    rvs(D) = mean(T.sunFPWMeanCutDiff(sIDays == uniqueNights(D)));
end
[weightsSI, bestCorrSI] = findWeightsCartesian(d, c, w, sI);
title('Sun:SI Weighting', 'Color', 'w', 'FontSize', 32)

tracker = table({'Sun'}, bestCorr, bestCorrSI, weightsSI, bestCorrSI, 'VariableNames', {'Target', 'vCon', 'ModelSI', 'WeightsSI', 'SI'});
save model.mat tracker weights dI dSigns cI cSigns wI wSigns

%Ranks lines by corr strength, tries averaging over 1:N, finds max N for
%corr with specified endpoint
function [b, series, I, signs] = optimizeComp(lineData, endpoint)
sdoComp = corr(lineData, endpoint);
[~, I] = sort(abs(sdoComp), 'descend');

I = I(~isnan(sdoComp(I)));

corrSorted = sdoComp(I);
I = I(abs(corrSorted) > .1);
corrSorted = corrSorted(abs(corrSorted) > .1);

signs = sign(corrSorted);

avgd = zeros(length(endpoint), length(I));

%Requiring at least 20 lines 
for i = 20:length(I)
    avgd(:, i) = mean(lineData(:, I(1:i)) .* signs(1:i)', 2);
end

[b, i] = max(corr(avgd, endpoint));
series = avgd(:, i);
I = I(1:i);
signs = signs(1:i);
end


