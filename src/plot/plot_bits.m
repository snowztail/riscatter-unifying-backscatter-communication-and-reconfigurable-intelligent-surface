clear; run('../setup'); run(strcat('config_', erase(mfilename, 'plot_'))); clc; close all;

%% * Load data
directory = strcat('../data/rate_', erase(mfilename, 'plot_'), '/');
data_load;

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
plotHandle(1) = plot(bits, rate(1, :) / log(2), 'DisplayName', 'Primary rate');
legend('Location', 'se'); grid on; box on; axis tight;
xticks(bits);
xlabel('Energy Discretization Bits');
ylabel('Primary Rate [bits/s/Hz]');
style_plot(plotHandle);
savefig('figures/rate_bits_primary.fig');
matlab2tikz('../../assets/simulation/rate_bits_primary.tex', 'extraaxisoptions', {'title style={font=\huge}', 'label style={font=\huge}', 'ticklabel style={font=\LARGE}', 'legend style={font=\LARGE}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=2}'});
% * Total backscatter rate
figure('Name', 'Average Total Backscatter Rate vs Energy Discretization Bits', 'Position', [0, 0, 500, 400]);
plotHandle(2) = plot(bits, rate(2, :) / log(2), 'DisplayName', 'Total backscatter rate');
legend('Location', 'se'); grid on; box on; axis tight;
xticks(bits);
xlabel('Energy Discretization Bits');
ylabel('Total Backscatter Rate [bits/BB]');
style_plot(plotHandle);
savefig('figures/rate_bits_backscatter.fig');
matlab2tikz('../../assets/simulation/rate_bits_backscatter.tex', 'extraaxisoptions', {'title style={font=\huge}', 'label style={font=\huge}', 'ticklabel style={font=\LARGE}', 'legend style={font=\LARGE}', 'scaled y ticks=false', 'y tick label style={/pgf/number format/.cd, fixed, precision=3}'});
