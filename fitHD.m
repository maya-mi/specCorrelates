load(strcat(starName, '/processed', starName, '.mat'))
load(strcat(starName, '/propError', starName, '.mat'))
load(strcat(starName, '/widthsOffsets', starName, '.mat'))
load(strcat(starName, '/times', starName, '.mat'))
ironA = ironA + offsets';

[nObs, nLines, nPix] = size(normOrders);
% finalDays = length(edges);
%Use custom fit for gaussian + vertical shift
mfit = @(b, x)(b(1) - b(2)*exp(-((x - b(3))/(sqrt(2)*b(4))).^2 * .5));

%Pre-designating fit vars
nNights = length(uniqueNights);
f = zeros(nNights, nLines, 4);
errFit = zeros(nNights, nLines, 4);
chisq = zeros(nNights, nLines);
reduced = zeros(nNights, nLines);


for i = 1:nNights
    for j = 1:nLines
        try
            thisNight = uniqueNights(i) == obsNights;
            fitX = squeeze(wavelengths(thisNight, j, :)) - ironA(j); 
            [~, lI] = min(abs(mean(fitX)  + 1.75 * widths(j) * sqrt(2)));
            [~, rI] = min(abs(mean(fitX)  - 1.75 * widths(j) * sqrt(2)));
            if isempty(lI)
                lI = 1;
            end
            if isempty(rI)
                rI = 1;
            end

            fitX = fitX(:, lI:rI);
            fitX = reshape(fitX, 1, numel(fitX));
            [fitX, I] = sort(fitX); %sorting to make bisectors easier
            fitY = squeeze(normOrders(thisNight, j, lI:rI)); 
            fitY = reshape(fitY, 1, numel(fitY));
            fitY = fitY(I);
            err = squeeze(propError(thisNight, j, lI:rI));
            err = reshape(err, 1, numel(err));
            errSq = err(I) .^2;

            [mn, loc] = min(fitY);
            b0 = [1, 1 - mn, fitX(loc), 0.035];
            [f(i, j, :), ~, ~, cov, ~, ~] = nlinfit(fitX, fitY, mfit, b0, 'Weights', (1 ./ errSq));
            errFit(i, j, :) = sqrt(diag(cov));
            model = mfit(f(i, j, :), fitX);
            chisq(i, j) = sum((fitY - model) .^2 ./ errSq);
            reduced(i, j) = chisq(i, j) ./ (length(fitX) - 4);
        catch
            f(i, j, :) = nan(1, 4);
            errFit(i, j, :) = nan(1, 4);
            chisq(i, j) = nan;
            reduced(i, j) = nan;
            
        end


    end
     
end

save(strcat(starName, '/fitResults', starName, '.mat'), 'f', 'errFit', 'chisq', 'reduced', 'ironA')
