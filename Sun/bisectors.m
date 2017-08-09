% load processedSun
% load propErrorSun
% load widthsOffsetsSun
% load fitResultsSun
% load ironA
% ironA = ironA + offsets';
load timesSun
[nObs, nLines, nPix] = size(normOrders);

nDays = length(uniqueNights);
%Bisector vars
nSteps = 20;
bisectY = zeros(nDays, nLines, nSteps * .9);
bisectX = zeros(nDays, nLines, nSteps * .9);
err = zeros(nDays, nLines, nSteps * .9);


for i = 1:nDays
            for j = 1:nLines

                %Coadding:
                    thisDay = uniqueNights(i) == obsNights;
                    fitX = squeeze(wavelengths(thisDay, j, :)) - ironA(j); 
                    [~, lI] = min(abs (mean(fitX) + 3 * sqrt(2) * widths(j)));
                    [~, rI] = min(abs (mean(fitX) - 3 * sqrt(2) * widths(j)));
                    %In case of nans: 
                    if isempty(lI)
                        lI = 1;
                    end
                    if isempty(rI)
                        rI = 31;
                    end
                  
                    fitX = fitX(:, lI:rI);
                    fitY = squeeze(normOrders(thisDay, j, lI:rI)); 
                    errX = squeeze(propError(thisDay, j, lI:rI));
                    ampStep = f(i, j, 2) / 10; %Don't care about the top 10%
                    [nO, nP] = size(fitX);
                    %Amplitude grid: * .9 to throw away top 10%:
                    ampGrid = linspace(f(i, j, 1) - f(i, j, 2), f(i, j, 1) - ampStep, nSteps * .9);
                    left = zeros(nO, length(ampGrid));
                    right = zeros(nO, length(ampGrid));
                    leftError = zeros(nO, length(ampGrid));
                    rightError = zeros(nO, length(ampGrid));

                    for k = 1:nO
                        [~, mid] = min(abs(fitX(k, :)));
                        left(k, :) = interp1(fitY(k, 1:mid), fitX(k, 1:mid), ampGrid);
                        right(k, :) = interp1(fitY(k, mid:end), fitX(k, mid:end), ampGrid);
                        leftError(k, :) = interp1(fitY(k, 1:mid), errX(k, 1:mid), ampGrid);
                        rightError(k, :) = interp1(fitY(k, mid:end), errX(k, mid:end), ampGrid);
                    end

                    lefts = nanmean(left);
                    rights = nanmean(right); 
                    bisectX(i, j, :) = mean([lefts; rights]);
                    bisectY(i, j, :) = ampGrid;
                    %Since y increments are constant, take dx/dy. Originally
                    %dividing by dy/dx --> multiply by dx/dy
                    lm = abs(gradient(lefts, mean(diff(ampGrid))))';
                    rm = abs(gradient(rights, mean(diff(ampGrid))))';
                    m = nanmean([lm, rm], 2)';
                    errF = mean(sqrt(leftError .^2 + rightError.^2));
                    err(i, j, :) = errF .* m / sqrt(2);
                
            end
end

bottom = round(.1*nSteps): round(.4*nSteps);
top = round(.6*nSteps):round(.9*nSteps);
bisectSpan = squeeze(nanmean(bisectX(:, :, top), 3)) - squeeze(nanmean(bisectX(:, :, bottom), 3));
spanErr = squeeze(nanmean(err(:, :, top), 3)) - squeeze(nanmean(err(:, :, bottom), 3));
%NOT SHOWING STRONG AS EXPECTED CORRELATION W/ CENTER. REVISIT. 
save('bisectSun.mat', 'bisectY', 'bisectX', 'bisectSpan', 'err', 'spanErr')

