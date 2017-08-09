load processedSun
load propErrorSun
load offsets
load widths
load fitResultsSun
load ironA
ironA = ironA + offsets;
[nObs, nLines, nPix] = size(normOrders);
fwhmGuess = 2*sqrt(2*log(2))*widths;

%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(-b(1)* (b(4)/2 * (b(4)/2) ./((x - b(3)).^2 + (b(4)/2)^2) + (1 - b(2))*exp(-((x - b(3))*2*sqrt(log(2))/b(4)).^2)) + b(5));
b0 = [.5 .5 0 .08 1];
bMin = [0 0 -.2 0 .5];
bMax = [1 1 .2 .3 1.5];
nParams = 6;
notEmpty = nums > 5;
%Pre-designating fit vars
pvF = zeros(sum(hasSDO), nLines, nParams);
errPVF = zeros(sum(hasSDO), nLines);

sdoI = 0;
index = 0;
for i = 1:max(D)
    if notEmpty(i) %Enforcing min number of daily exposures
        sdoI = sdoI + 1;
        if hasSDO(sdoI)
            index = index + 1;
            for j = 1:nLines
                
                    fitX = squeeze(wavelengths(i == D, j, :)) - ironA(j); 
                    [~, ll] = min(abs(mean(fitX) +  3.5 * widths(j) *sqrt(2)));
                    [~, rr] = min(abs(mean(fitX) -3.5 * widths(j) *sqrt(2)));
                    
                    if isempty(ll)
                        ll = 1;
                    end
                    if isempty(rr)
                        rr = 1;
                    end
                    
                    fitX = fitX(:, ll:rr);
                    fitX = reshape(fitX, 1, numel(fitX));
                    [fitX, I] = sort(fitX); 
                     
                    fitY = squeeze(normOrders(i == D, j, ll:rr)); 
                    fitY = reshape(fitY, 1, numel(fitY));
                    fitY = fitY(I);

                    [pvF(index, j, :), errPVF] = lsqcurvefit(mfit, b0, fitX, fitY,bMin, bMax); 
                
            end
        end
               
    end
end

save('pvfNew.mat', 'pvF', 'errPVF')

