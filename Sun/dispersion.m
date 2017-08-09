load('processedSun.mat')
load('propErrorSun.mat')
load('offsets.mat')
load widths.mat
ironA = ironA + offsets;
[nObs, nLines, nPix] = size(normOrders);
finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(-1 / pi * b(1) ./ ((x - b(2)).^2 + b(1)^2)*b(4) + b(3));
b0 = [.001, 0, 1, .5];
nParams = 4;
notEmpty = nums > 5;
%Pre-designating fit vars
fL = zeros(sum(hasSDO), nLines, nParams);
fR = zeros(sum(hasSDO), nLines, nParams);
fB = zeros(sum(hasSDO), nLines, nParams);
errFitL = zeros(sum(hasSDO), nLines, nParams);
errFitR = zeros(sum(hasSDO), nLines, nParams);
errFitB = zeros(sum(hasSDO), nLines, nParams);
reducedL = zeros(sum(hasSDO), nLines);
reducedR = zeros(sum(hasSDO), nLines);
reducedB = zeros(sum(hasSDO), nLines);
sdoI = 0;
index = 0;
for i = 1:max(D)
    if notEmpty(i) %Enforcing min number of daily exposures
        sdoI = sdoI + 1;
        if hasSDO(sdoI)
            index = index + 1;
            for j = 1:nLines
                try 
                    fitX = squeeze(wavelengths(i == D, j, :)) - ironA(j);
                    [~, lI] = min(abs(mean(fitX) + 1.75 * widths(j)*sqrt(2)));
                    [~, rI] = min(abs(mean(fitX) - 1.75 * widths(j)*sqrt(2)));
                    [~, ll] = min(abs(mean(fitX) + 3.5 * widths(j)*sqrt(2)));
                    [~, rr] = min(abs(mean(fitX) - 3.5 * widths(j)*sqrt(2)));

                    if isempty(ll)
                        ll = 1;
                    end
                    if isempty(rr)
                        rr = 31;
                    end
                
                    fitXL = fitX(:, ll:lI);
                    fitXR = fitX(:, rI:rr);
                    fitXL = reshape(fitXL, 1, numel(fitXL));
                    fitXR = reshape(fitXR, 1, numel(fitXR));
                    [fitXL, IL] = sort(fitXL); 
                    [fitXR, IR] = sort(fitXR); 
                    fitYL = squeeze(normOrders(i == D, j, ll:lI)); 
                    fitYL = reshape(fitYL, 1, numel(fitYL));
                    fitYL = fitYL(IL);
                    fitYR = squeeze(normOrders(i == D, j, rI:rr)); 
                    fitYR = reshape(fitYR, 1, numel(fitYR));
                    fitYR = fitYR(IR);
                    
                    errL = squeeze(propError(i == D, j, ll:lI));
                    errL = reshape(errL, 1, numel(errL));
                    errLSq = errL(IL) .^2;
                    errR = squeeze(propError(i == D, j, rI:rr));
                    errR = reshape(errR, 1, numel(errR));
                    errRSq = errR(IR) .^2;
                    errBSq = [errLSq, errRSq];
                   
                    [fL(index, j, :), ~, ~, covL, ~, ~] = nlinfit(fitXL, fitYL, mfit, b0, 'Weights', (1 ./ errLSq));
                    [fR(index, j, :), ~, ~, covR, ~, ~] = nlinfit(fitXR, fitYR, mfit, b0, 'Weights', (1 ./ errRSq));
                    [fB(index, j, :), ~, ~, covB, ~, ~] = nlinfit([fitXL, fitXR], [fitYL, fitYR], mfit, b0, 'Weights', (1 ./ errBSq));
                    errFitL(index, j, :) = sqrt(diag(covL));
                    errFitR(index, j, :) = sqrt(diag(covR));
                    errFitB(index, j, :) = sqrt(diag(covB));
                    modelL = mfit(fL(index, j, :), fitXL);
                    chisqL = sum((fitYL - modelL) .^2 ./ errLSq);
                    reducedL(index, j) = chisqL ./ (length(fitXL) - nParams);
                    modelR = mfit(fR(index, j, :), fitXR);
                    chisqR = sum((fitYR - modelR) .^2 ./ errRSq);
                    reducedR(index, j) = chisqR ./ (length(fitXR) - nParams);
                    modelB = mfit(fB(index, j, :), [fitXL, fitXR]);
                    chisqB = sum(([fitYL, fitYR] - modelB) .^2 ./ errBSq);
                    reducedB(index, j) = chisqB ./ (length(errBSq) - nParams);
                catch
                end
            end
        end
               
    end
end

save('dispersionSun.mat', 'fL', 'fR', 'fB', 'reducedL', 'reducedR', 'reducedB', 'errFitL', 'errFitR', 'errFitB', 'ironA')

