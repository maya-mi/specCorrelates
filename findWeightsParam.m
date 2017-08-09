function[weightsBest, best] = findWeightsParam(d, c, w, endpoint)
theta = 0:.05:2*pi;
phi = 0:.05:pi/2;

corrs = zeros(length(theta), length(phi));
best = 0;
weightsBest = [0 0 0];

for i = 1:length(theta)
    for j = 1:length(phi)
        A = cos(theta(i))*sin(phi(j));
        B = sin(theta(i))*sin(phi(j));
        C = cos(phi(j));
        combo = A*d + B*c + C*w;
        corrs(i, j) = corr(combo, endpoint);
        if abs(corrs(i, j)) > abs(best)
            best = corrs(i, j);
            weightsBest = [A B C];
        end
    end
end

figure; imagesc(corrs)
xlabel('Phi')
ylabel('Theta')
axis equal
box off
ax = gca;
set(gcf,'DefaultTextColor','w');

ax.FontSize = 18;
ax.XColor = 'w';
ax.YColor = 'w';
c = colorbar;
c.Color = 'w';
axis equal

set(gca, 'Color', 'none')
xlim([0, length(phi)])
ylim([0, length(theta)])
xTicks = 3;
yTicks = 6;
phiLabel = cellfun(@num2str, num2cell(linspace(min(phi), max(phi), xTicks)), 'UniformOutput', 0);
thetaLabel = cellfun(@num2str, num2cell(linspace(min(theta), max(theta), yTicks)), 'UniformOutput', 0);
set(gca, 'XTick', linspace(1, length(phi), xTicks))
set(gca, 'YTick', linspace(1, length(theta), yTicks))
set(gca,'YTickLabel', thetaLabel)
set(gca,'XTickLabel', phiLabel)
end