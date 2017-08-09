load(strcat(starName, '/processed', starName, '.mat'))
load(strcat(starName, '/propError', starName, '.mat'))
load(strcat(starName, '/times', starName, '.mat'))
load widths.mat

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


for i = 1:length(uniqueNights)
    for j = 1:nLines
        try
            thisNight = uniqueNights(i) == obsNights;
            fitX = squeeze(wavelengths(thisNight, j, :)) - ironA(j); 
            [~, lI] = min(abs(mean(fitX)  + 2 * widths(j) * sqrt(2)));
            [~, rI] = min(abs(mean(fitX)  - 2 * widths(j) * sqrt(2)));
            if isempty(lI)
                lI = 1;
            end
            if isempty(rI)
                rI = 31;
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

widths = nanmean(abs(f(:, :, 4)));
widths(widths > .1) = .1;

offsets = nanmean(f(:, :, 3));
offsets(offsets > .1) = .1;
offsets(offsets < -.1) = -.1;
save(strcat(starName, '/widthsOffsets', starName, '.mat'), 'widths', 'offsets');
%save('fitResultsHD.mat', 'f', 'errFit', 'chisq', 'reduced', 'ironA')

