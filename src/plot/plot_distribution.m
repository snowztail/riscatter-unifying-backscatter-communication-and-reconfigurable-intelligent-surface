clear; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
directory = strcat('../data/region_', erase(mfilename, 'plot_'), '/');
data_load;

%% * Average over instances and retrieve rate regions
region = cell(nVariables, 1);
for iVariable = 1 : nVariables
	rate = zeros(2, nWeights + 3);
	for iWeight = 1 : nWeights
		rate(:, iWeight) = mean(cat(2, Result(iVariable, iWeight, :).rate), 2);
	end
	[rate(1, nWeights + 1), rate(2, nWeights + 2)] = deal(max(rate(1, :)), max(rate(2, :)));
	region{iVariable} = rate(:, convhull(transpose(rate)));
end
save(strcat('../data/region_', erase(mfilename, 'plot_')));

%% * Draw primary-(sum-)backscatter rate regions
figureHandle = figure('Name', 'Average Primary-(Sum-)Backscatter Rate Region by Different Input Distribution Design', 'Position', [0, 0, 500, 400]);
plotHandle = gobjects(nVariables, 1);
hold all;
for iVariable = 1 : nVariables
	plotHandle(iVariable) = plot(region{iVariable}(1, :) / log(2), region{iVariable}(2, :) / log(2));
end
hold off; legend({'Cooperation', 'Exhaustion', 'KKT', 'Equiprobable', 'Marginalization', 'Decomposition', 'Randomization'}, 'Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Backscatter Rate [bits/BB]');
xlim([6, Inf]);
style_plot(plotHandle);
magnifyOnFigure(figureHandle, 'Mode', 'interactive', 'initialPositionMagnifier', [330 205 10 10], 'initialPositionSecondaryAxes', [150 250 100 100], 'secondaryAxesFaceColor', [0.9 0.9 0.9], 'secondaryAxesXLim', [6.33 6.37], 'secondaryAxesYLim', [0.079 0.084]);
savefig(strcat('figures/region_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/simulation/region_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', {'title style={font=\LARGE}', 'label style={font=\LARGE}', 'ticklabel style={font=\Large}', 'legend style={font=\Large}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=2}'});
