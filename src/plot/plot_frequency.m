clear; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
directory = strcat('../data/region_', erase(mfilename, 'plot_'), '/');
data_load;

%% * Average over instances and retrieve rate regions
region = cell(nVariables, 1);
for iVariable = 1 : nVariables
	rate = zeros(2, nWeights + 3);
	for iWeight = 1 : nWeights
		rate(:, iWeight) = mean(cat(2, Result(iVariable, iWeight, :).rate), 2, 'omitnan');
	end
	[rate(1, nWeights + 1), rate(2, nWeights + 2)] = deal(max(rate(1, :)), max(rate(2, :)));
	region{iVariable} = rate(:, convhull(transpose(rate)));
end
save(strcat('../data/region_', erase(mfilename, 'plot_')));

%% * Draw primary-(sum-)backscatter rate regions
figure('Name', 'Average Primary-(Sum-)Backscatter Rate Region vs Carrier Frequency', 'Position', [0, 0, 500, 400]);
object = gobjects(nVariables, 1);
hold all;
for iVariable = 1 : nVariables
	frequency = Variable(iVariable).frequency;
	object(iVariable) = plot(region{iVariable}(1, :), region{iVariable}(2, :), 'DisplayName', strcat('$f = ', num2str(frequency / 1e6), '$ MHz'));
end
hold off; legend('Location', 'se'); grid on; box on; axis tight;
set(gca, 'YScale', 'log');
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Total Backscatte Rate [$\mu$ bits/backscatter symbol]');
ylim([1e-7, 1e-1]);
style_plot(object);
savefig(strcat('figures/region_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/region_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
