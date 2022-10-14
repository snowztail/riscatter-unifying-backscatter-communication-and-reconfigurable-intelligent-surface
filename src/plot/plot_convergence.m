clear; run('../setup'); load('../data/wsr_convergence'); clc; close all;

%% * Draw convergence curves
figure('Name', 'Convergence behavior of the proposed algorithms', 'Position', [0, 0, 500, 400]);
convergencePlot = tiledlayout(3, 1, 'tilespacing', 'compact');
plotHandle = gobjects(1, 3);

% * KKT (first call)
nexttile;
plotHandle(1) = plot(0 : 1e1 : 1.2e2, convergence.kkt(1 : 1e1 : 1.2e2 + 1), 'DisplayName', 'KKT (First Call)');
hold off; legend('Location', 'se'); grid on; box on;

% * PGD (first call)
nexttile;
plotHandle(2) = plot(0 : length(convergence.pgd) - 1, convergence.pgd, 'DisplayName', 'PGD (First Call)');
hold off; legend('Location', 'se'); grid on; box on;
ylabel('Weighed Sum-Rate');

% * BCD
nexttile;
plotHandle(3) = plot(0 : length(convergence.bcd) - 1, convergence.bcd, 'DisplayName', 'BCD');
hold off; legend('Location', 'se'); grid on; box on;
xlabel('Number of Iterations');
style_plot(plotHandle);

savefig('figures/wsr_convergence');
matlab2tikz('../../assets/simulation/wsr_convergence.tex', 'extraaxisoptions', ['title style={font=\Large}, ' 'label style={font=\Large}, ' 'ticklabel style={font=\large}, ' 'legend style={font=\large}, ' 'yticklabel=\pgfkeys{/pgf/number format/.cd,fixed,precision=3}\pgfmathprintnumber{\tick}']);
