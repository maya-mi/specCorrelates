load starInfo.mat
load model.mat

for starModelCounter = 1:length(stars)
    starName = stars{starModelCounter};
    if ~strcmp(starName, 'Sun')
        load(strcat(starName, '/fitResults', starName, '.mat'))
        load(strcat(starName, '/times', starName, '.mat'))
        nNights = length(uniqueNights);
        %Unable to recognize .rdb format? Place desired stars info into
        %starsInfo folder, and convert manually to .txt (for now)
        starData = readtable(strcat('starsInfo/', starName, '_harpn.txt'), 'HeaderLines', 1);
        refTimes = datetime(2.4e+6 + starData{:, 1}, 'ConvertFrom', 'juliandate');
        refNights = dateshift(refTimes - hours(12), 'start', 'day');
        rvs = starData{:, 2};
        sI = starData{:, 8};
        avgSI = zeros(nNights, 1);
        avgRVs = zeros(nNights, 1);
        for nCounter = 1:nNights
            avgSI(nCounter) = mean(sI(uniqueNights(nCounter) == refNights));
            avgRVs(nCounter) = mean(rvs(uniqueNights(nCounter) == refNights));
        end
    
        wData = squeeze(f(:, wI, 4));
        wData = abs(wData);
        errSqW = squeeze(errFit(:, wI, 4)).^2;
        errSqW(wData > .1 | errSqW == inf) = nan;
        wErrWeights = normWeights(errSqW, 100);
        widthCombo = nansum(wData.* wErrWeights .* wSigns', 2);
        widthCombo(widthCombo == 0) = nan;

        cData = squeeze(f(:, cI, 3));
        errSqC = squeeze(errFit(:, cI, 3)).^2;
        errSqC( errSqC == 0 | errSqC == inf | cData > .05) = nan;
        cErrWeights = normWeights(errSqC, 100);
        centerCombo = nansum(cData .* cErrWeights .* cSigns', 2);
        centerCombo(centerCombo == 0) = nan;
        
        dData = squeeze(f(:, dI, 2) ./ f(:, dI, 1));
        errSqD = squeeze(errFit(:, dI, 1)).^2 + squeeze(errFit(:, dI, 2)).^2;
        errSqD(errSqD == 0 | errSqD == inf) = nan;
        dErrWeights = normWeights(errSqD, 100);
        depthCombo = nansum(dData .* dErrWeights .* dSigns', 2);
        depthCombo(depthCombo == 0) = nan;

        [d, dsc, ~] = normSeries(depthCombo);
        [w, wsc, ~] = normSeries(widthCombo);
        [c, csc, ~] = normSeries(centerCombo);

        model = weights(1)*d + weights(2)*c + weights(3)*w;
        nanfilter = ~isnan(model) & ~isnan(avgSI);
        modelCorr = corr(model(nanfilter), avgSI(nanfilter));
        
        
        [weightsSI, bestCorrSI] = findWeightsCartesian(d(nanfilter), c(nanfilter), w(nanfilter), avgSI(nanfilter));
        title(strcat(starName, ':SI corr'))
        tracker{starModelCounter, 1} = {starName}; 
        tracker{starModelCounter, 2} = nan;
        tracker{starModelCounter, 3} = modelCorr;
        tracker{starModelCounter, 4} = weightsSI;
        tracker{starModelCounter, 5} = bestCorrSI;
    end


end

save testedModel.mat tracker


function [nW] = normWeights(invErrSq, maxReps)
    
    wSum = nansum(invErrSq, 2);
    wSum(wSum == inf) = nan;
    nW = invErrSq ./ wSum;
    if maxReps > 0
        if ~isempty(find(nW > .4, 1))
            nW(nW > .2) = .35;
            nW = normWeights(nW, maxReps - 1);
        end
    else 
    end
end
    