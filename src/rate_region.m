clear; cvx_clear; clc; close all; setup;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
nWeights = 3e1;
weightSet = [linspace(0, 0.25, nWeights - 1), 1];
reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;
nBins = 2 ^ 8;
constellation = qammod(0 : nStates - 1, nStates);
confidenceScore = 10;
tolerance = 1e-6;
precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Compute expected received power per primary symbol
indexCombination = index_combination(nTags, nStates);
inputCombination = transpose(constellation(indexCombination));
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(reflectRatio) * inputCombination(iInput, :) * cascadedChannel;
end
signalPower = abs(equivalentChannel * precoder) .^ 2;
snr = signalPower / noisePower;
receivedPower = signalPower + noisePower;

% * Obtain threshold candidates within empirical interval based on Chebyshev's inequality
lowerBound = symbolRatio * min(receivedPower) - confidenceScore * sqrt(symbolRatio * min(receivedPower));
upperBound = symbolRatio * max(receivedPower) + confidenceScore * sqrt(symbolRatio * max(receivedPower));
if lowerBound <= 0
	thresholdCandidate = [linspace(0, upperBound, nBins), inf];
else
	thresholdCandidate = [0, linspace(lowerBound, upperBound, nBins - 1), inf];
end

% * Discretize continuous output and remaps to DMC
dmc = channel_discretization(thresholdCandidate, receivedPower, symbolRatio);

% * Evaluate rate regions
rateCooperationSmawk = zeros(nWeights, 2);
rateCooperationBisection = zeros(nWeights, 2);
rateCooperationMl = zeros(nWeights, 2);

rateExhaustionSmawk = zeros(nWeights, 2);
rateExhaustionBisection = zeros(nWeights, 2);
rateExhaustionMl = zeros(nWeights, 2);

rateKktSmawk = zeros(nWeights, 2);
rateKktBisection = zeros(nWeights, 2);
rateKktMl = zeros(nWeights, 2);

rateMarginalizationSmawk = zeros(nWeights, 2);
rateMarginalizationBisection = zeros(nWeights, 2);
rateMarginalizationMl = zeros(nWeights, 2);

for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	rateCooperationSmawk(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'cooperation', 'Threshold', 'smawk');
	rateCooperationBisection(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'cooperation', 'Threshold', 'bisection');
	rateCooperationMl(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'cooperation', 'Threshold', 'ml');

	rateExhaustionSmawk(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'exhaustion', 'Threshold', 'smawk');
	rateExhaustionBisection(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'exhaustion', 'Threshold', 'bisection');
	rateExhaustionMl(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'exhaustion', 'Threshold', 'ml');

	rateKktSmawk(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'kkt', 'Threshold', 'smawk');
	rateKktBisection(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'kkt', 'Threshold', 'bisection');
	rateKktMl(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'kkt', 'Threshold', 'ml');

	rateMarginalizationSmawk(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'marginalization', 'Threshold', 'smawk');
	rateMarginalizationBisection(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'marginalization', 'Threshold', 'bisection');
	rateMarginalizationMl(iWeight, :) = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, 'Input', 'marginalization', 'Threshold', 'ml');
end

save('data/rate_region.mat');

figureHandle = figure('name', 'Achievable rate regions by different input and threshold design', 'position', [0, 0, 500, 400]);
rateCooperationSmawk(end + 2, 2) = max(rateCooperationSmawk(:, 2));
rateCooperationBisection(end + 2, 2) = max(rateCooperationBisection(:, 2));
rateCooperationMl(end + 2, 2) = max(rateCooperationMl(:, 2));

rateExhaustionSmawk(end + 2, 2) = max(rateExhaustionSmawk(:, 2));
rateExhaustionBisection(end + 2, 2) = max(rateExhaustionBisection(:, 2));
rateExhaustionMl(end + 2, 2) = max(rateExhaustionMl(:, 2));

rateKktSmawk(end + 2, 2) = max(rateKktSmawk(:, 2));
rateKktBisection(end + 2, 2) = max(rateKktBisection(:, 2));
rateKktMl(end + 2, 2) = max(rateKktMl(:, 2));

rateMarginalizationSmawk(end + 2, 2) = max(rateMarginalizationSmawk(:, 2));
rateMarginalizationBisection(end + 2, 2) = max(rateMarginalizationBisection(:, 2));
rateMarginalizationMl(end + 2, 2) = max(rateMarginalizationMl(:, 2));

plotHandle = gobjects(3, 4);
legendHandle = cell(3, 4);
hold all;

plotHandle(1, 1) = plot(rateCooperationSmawk(convhull(rateCooperationSmawk), 1), rateCooperationSmawk(convhull(rateCooperationSmawk), 2));
plotHandle(2, 1) = plot(rateCooperationBisection(convhull(rateCooperationBisection), 1), rateCooperationBisection(convhull(rateCooperationBisection), 2));
plotHandle(3, 1) = plot(rateCooperationMl(convhull(rateCooperationMl), 1), rateCooperationMl(convhull(rateCooperationMl), 2));
legendHandle{1, 1} = 'Cooperation-SMAWK';
legendHandle{2, 1} = 'Cooperation-Bisection';
legendHandle{3, 1} = 'Cooperation-ML';

plotHandle(1, 2) = plot(rateExhaustionSmawk(convhull(rateExhaustionSmawk), 1), rateExhaustionSmawk(convhull(rateExhaustionSmawk), 2));
plotHandle(2, 2) = plot(rateExhaustionBisection(convhull(rateExhaustionBisection), 1), rateExhaustionBisection(convhull(rateExhaustionBisection), 2));
plotHandle(3, 2) = plot(rateExhaustionMl(convhull(rateExhaustionMl), 1), rateExhaustionMl(convhull(rateExhaustionMl), 2));
legendHandle{1, 2} = 'Exhaustion-SMAWK';
legendHandle{2, 2} = 'Exhaustion-Bisection';
legendHandle{3, 2} = 'Exhaustion-ML';

plotHandle(1, 3) = plot(rateKktSmawk(convhull(rateKktSmawk), 1), rateKktSmawk(convhull(rateKktSmawk), 2));
plotHandle(2, 3) = plot(rateKktBisection(convhull(rateKktBisection), 1), rateKktBisection(convhull(rateKktBisection), 2));
plotHandle(3, 3) = plot(rateKktMl(convhull(rateKktMl), 1), rateKktMl(convhull(rateKktMl), 2));
legendHandle{1, 3} = 'KKT-SMAWK';
legendHandle{2, 3} = 'KKT-Bisection';
legendHandle{3, 3} = 'KKT-ML';

plotHandle(1, 4) = plot(rateMarginalizationSmawk(convhull(rateMarginalizationSmawk), 1), rateMarginalizationSmawk(convhull(rateMarginalizationSmawk), 2));
plotHandle(2, 4) = plot(rateMarginalizationBisection(convhull(rateMarginalizationBisection), 1), rateMarginalizationBisection(convhull(rateMarginalizationBisection), 2));
plotHandle(3, 4) = plot(rateMarginalizationMl(convhull(rateMarginalizationMl), 1), rateMarginalizationMl(convhull(rateMarginalizationMl), 2));
legendHandle{1, 4} = 'Marginalization-SMAWK';
legendHandle{2, 4} = 'Marginalization-Bisection';
legendHandle{3, 4} = 'Marginalization-ML';

hold off;
grid minor;
box on;
legend(legendHandle(:), 'Location', 'southwest');
xlabel('Primary rate [nats/s/Hz]');
ylabel('Baskcatter sum rate [nats/channel use]');
plot_style_group(plotHandle(:), 3);
magnifyOnFigure(figureHandle, 'Mode', 'interactive', 'initialPositionMagnifier', [360 250 20 20], 'initialPositionSecondaryAxes', [300 65 100 100], 'secondaryAxesFaceColor', [0.9 0.9 0.9]);
