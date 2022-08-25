clear; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
directory = strcat('../data/rate_', erase(mfilename, 'plot_'), '/');
data_load;
% Result = load(strcat(directory, 'instance_0'), 'Result').Result;

%% * Obtain average total backscatter rate
bits = cat(2, Variable.nBits);
rate = zeros(2, nVariables);
for iVariable = 1 : nVariables
	rate(:, iVariable) = mean(cat(2, Result(iVariable, 1, :).rate), 2);
end
save(strcat('../data/rate_', erase(mfilename, 'plot_')));

%% * Draw rate curves
% * Primary rate
figure('Name', 'Average Primary Rate vs Energy Discretization Bits', 'Position', [0, 0, 500, 400]);
object(1) = plot(bits, rate(1, :) / log(2), 'DisplayName', 'Primary rate');
legend('Location', 'se'); grid on; box on; axis tight;
xticks(bits);
% ylim([6, 10]);
xlabel('Energy Discretization Bits');
ylabel('Primary Rate [bits/s/Hz]');
style_plot(object);
savefig('figures/rate_bits_primary.fig');
matlab2tikz('../../assets/simulation/rate_bits_primary.tex', 'extraaxisoptions', {'title style={font=\huge}', 'label style={font=\huge}', 'ticklabel style={font=\LARGE}', 'legend style={font=\LARGE}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=2}'});
% * Total backscatter rate
figure('Name', 'Average Total Backscatter Rate vs Energy Discretization Bits', 'Position', [0, 0, 500, 400]);
object(2) = plot(bits, rate(2, :) / log(2), 'DisplayName', 'Total backscatter rate');
legend('Location', 'se'); grid on; box on; axis tight;
xticks(bits);
% ylim([0, 1]);
xlabel('Energy Discretization Bits');
ylabel('Total Backscatte Rate [bits/BSP]');
style_plot(object);
savefig('figures/rate_bits_backscatter.fig');
matlab2tikz('../../assets/simulation/rate_bits_backscatter.tex', 'extraaxisoptions', {'title style={font=\huge}', 'label style={font=\huge}', 'ticklabel style={font=\LARGE}', 'legend style={font=\LARGE}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=3}'});
