clear; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
distribution = horzcat(load(strcat('../data/distribution_', erase(mfilename, 'plot_')), 'Result').Result.distribution);

%% * Draw reflection state distributions
figure('Name', 'Tag Input Distribution vs Weight', 'Position', [0, 0, 500, 400]);
object = gobjects(nWeights, 1);
hold all;
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	object(iWeight) = plot(distribution(:, iWeight), 'DisplayName', strcat('$\rho = ', num2str(weight), '$'));
end
hold off; legend('Location', 'ne'); grid on; box on;
xlabel('Reflection State');
ylabel('Probability Distribution');
xticks(1 : nStates);
yticks(0 : 0.2 : 1);
style_plot(object);
savefig(strcat('figures/distribution_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/simulation/distribution_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
