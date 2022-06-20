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
figure('Name', 'Average Primary-(Sum-)Backscatter Rate Region by Different Decision Threshold Design', 'Position', [0, 0, 500, 400]);
object = gobjects(nVariables, 1);
hold all;
for iVariable = 1 : nVariables
	object(iVariable) = plot(region{iVariable}(1, :), 1e3 * region{iVariable}(2, :));
end
hold off; legend({'SMAWK', 'DP', 'Bisection', 'ML'}, 'Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [nats/s/Hz]');
ylabel('(Sum-)Backscatte Rate [$\mu$ nats/backscatter symbol duration]');
xlim([1.3, inf]);
style_plot(object);
savefig(strcat('figures/region_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/region_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);