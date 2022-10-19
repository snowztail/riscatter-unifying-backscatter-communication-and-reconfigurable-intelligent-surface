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
plotHandle(5) = plot(regionRiscatter(1, :) / log(2), regionRiscatter(2, :) / log(2), 'Color', '#77AC30', 'Marker', '^', 'DisplayName', 'RIScatter');
plotHandle(4) = scatter(rateRis(1, :) / log(2), rateRis(2, :) / log(2), 200, [0.4940 0.1840 0.5560], 'Marker', 'x', 'DisplayName', 'RIS');
plotHandle(3) = scatter(rateSr(1, :) / log(2), rateSr(2, :) / log(2), 200, [0.9290 0.6940 0.1250], 'Marker', 's', 'DisplayName', 'SR');
plotHandle(2) = scatter(rateAmbc(1, :) / log(2), rateAmbc(2, :) / log(2), 200, [0.8500 0.3250 0.0980], 'Marker', '+', 'DisplayName', 'AmBC');
plotHandle(1) = scatter(rateBbc(1, :) / log(2), rateBbc(2, :) / log(2), 200, [0 0.4470 0.7410], 'Marker', 'o', 'DisplayName', 'BBC');
hold off; legend(plotHandle, 'Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Total Backscatter Rate [bits/BSP]');

savefig('figures/region_comparison');
matlab2tikz('../../assets/simulation/region_comparison.tex', 'extraaxisoptions', {'title style={font=\LARGE}', 'label style={font=\LARGE}', 'ticklabel style={font=\Large}', 'legend style={font=\Large}', 'reverse legend', 'every axis plot/.append style={line width=2pt}'});
