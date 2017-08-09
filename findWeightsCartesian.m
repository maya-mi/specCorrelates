function[weightsBest, best] = findWeightsCartesian(d, c, w, endpoint)
nSteps = 100;
x = linspace(-1, 1, nSteps);
xs = repmat(x, nSteps, 1);
corrs = zeros(nSteps);
best = 0;

weightsBest = [0 0 0];
zs = zeros(nSteps);

for i = 1: nSteps
    for j = 1:nSteps
        if x(i) ^2 + x(j)^2 <= 1
            zs(i, j) = sqrt(1 - x(i)^2 - x(j)^2);
            combo = d*x(i) + c*x(j) + w*zs(i, j);
            corrs(i, j) = corr(combo, endpoint);
            if abs(corrs(i, j)) > abs(best)
                best = corrs(i, j);
                weightsBest = [x(i), x(j), zs(i, j)];
            end
        else
            xs(i, j) = nan;
            zs(i, j) = nan;
        end
    end
end


figure; surf(xs, xs', zs, corrs);
ax = gca;
set(gcf,'DefaultTextColor','w');

ax.FontSize = 24;
ax.XColor = 'w';
ax.YColor = 'w';
xlabel('Depth')
ylabel('Center')
zlabel('Width')
c = colorbar;
c.Color = 'w';
axis equal
xlim([-1 1])
ylim([-1, 1])

set(gca, 'Color', 'none')
box off
view(2)

end