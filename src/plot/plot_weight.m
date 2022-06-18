clear; clear block_coordinate_descent beamforming_pgd; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
distribution = permute(load(strcat('../data/distribution_', erase(mfilename, 'plot_')), 'Result').Result.distribution, [1 3 2]);

%% * Draw reflection state distributions
figure('Name', 'Reflection State Distribution vs Weight', 'Position', [0, 0, 500, 400]);
object = gobjects(nWeights, 1);
hold all;
for iWeight = 1 : nWeights
	object(iWeight) = plot(distribution(:, iWeight), 'DisplayName', strcat('$\rho = ', num2str(weightSet(iWeight)), '$'));
end
hold off; legend; grid minor; box on;
xlabel('Reflection State');
ylabel('Probability Distribution');
xticks(1 : nStates);
style_plot(object);
savefig(strcat('figures/distribution_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/distribution_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
