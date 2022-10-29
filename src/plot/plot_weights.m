clear; run('../setup'); load('distribution_weights'); clc; close all;
distribution = horzcat(Result.distribution);

%% * Draw reflection state distributions
figure('Name', 'Tag Input Distribution vs Weight', 'Position', [0, 0, 500, 200]);
plotHandle = gobjects(nWeights, 1);
hold all;
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	plotHandle(iWeight) = plot(distribution(:, iWeight), 'DisplayName', strcat('$\rho = ', num2str(weight), '$'));
end
hold off; legend('Location', 'ne'); grid on; box on;
xlabel('Reflection State');
ylabel('Probability\\Distribution');
xticks(1 : nStates);
yticks(0 : 0.2 : 1);
style_plot(plotHandle);
savefig(strcat('figures/distribution_', erase(mfilename, 'plot_')));
matlab2tikz(strcat('../../assets/simulation/distribution_', erase(mfilename, 'plot_'), '.tex'), 'extraaxisoptions', {'align=center', 'title style={font=\LARGE}', 'label style={font=\LARGE}', 'ticklabel style={font=\Large}', 'legend style={font=\Large}'});
