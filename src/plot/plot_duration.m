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
figure('Name', 'Average Primary-(Sum-)Backscatter Rate Region vs Backscatter/Primary Symbol Duration Ratio', 'Position', [0, 0, 500, 400]);
object = gobjects(nVariables, 1);
hold all;
for iVariable = 1 : nVariables
	symbolRatio = Variable(iVariable).symbolRatio;
	object(iVariable) = plot(region{iVariable}(1, :) / log(2), region{iVariable}(2, :) / log(2) / symbolRatio, 'DisplayName', strcat('$N = ', num2str(symbolRatio), '$'));
end
hold off; legend('Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Total Backscatte Rate [bits/PSP]');
xlim([9, inf]);
style_plot(object);
savefig(strcat('figures/region_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/simulation/region_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', {'title style={font=\huge}', 'label style={font=\huge}', 'ticklabel style={font=\LARGE}', 'legend style={font=\LARGE}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=3}'});
