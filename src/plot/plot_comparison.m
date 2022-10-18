clear; run('../setup'); load('../data/region_comparison'); clc; close all;

%% * Retrieve RIScatter rate region
rate = zeros(2, nWeights + 3);
for iWeight = 1 : nWeights
	rate(:, iWeight) = ResultRiscatter(iWeight).rate;
end
[rate(1, nWeights + 1), rate(2, nWeights + 2)] = deal(max(rate(1, :)), max(rate(2, :)));
regionRiscatter = rate(:, convhull(transpose(rate)));

%% * Plot achievable rates
figure('Name', 'Achievable Rates of Different Scattering Applications', 'Position', [0, 0, 500, 400]);
plotHandle = gobjects(1, 5);
hold all;
plotHandle(1) = scatter(rateBbc(1, :) / log(2), rateBbc(2, :) / log(2), 'Marker', 'o', 'DisplayName', 'BBC');
plotHandle(2) = scatter(rateAmbc(1, :) / log(2), rateAmbc(2, :) / log(2), 'Marker', '+', 'DisplayName', 'AmBC');
plotHandle(3) = scatter(rateSr(1, :) / log(2), rateSr(2, :) / log(2), 'Marker', 'p', 'DisplayName', 'SR');
plotHandle(4) = scatter(rateRis(1, :) / log(2), rateRis(2, :) / log(2), 'Marker', 'x', 'DisplayName', 'RIS');
plotHandle(5) = plot(regionRiscatter(1, :) / log(2), regionRiscatter(2, :) / log(2), 'Marker', '^', 'DisplayName', 'RIScatter');
uistack(plotHandle(4), 'top');
hold off; legend(plotHandle, 'Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Total Backscatter Rate [bits/BSP]');

savefig('figures/region_comparison');
matlab2tikz('../../assets/simulation/region_comparison.tex', 'extraaxisoptions', ['title style={font=\Large}, ' 'label style={font=\Large}, ' 'ticklabel style={font=\large}, ' 'legend style={font=\large}']);

