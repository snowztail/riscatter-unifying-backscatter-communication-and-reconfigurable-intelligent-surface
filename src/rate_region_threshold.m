setup; clear; cvx_clear; clc; close all;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
nWeights = 2e1;
weightSet = [linspace(0, 0.1, nWeights - 1), 1];
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

% * Initialize input distribution, detection threshold, and DMTC
inputDistribution = ones(nTags, nStates) / nStates;
combinationDistribution = combination_distribution(inputDistribution);
equivalentDistribution = prod(combinationDistribution, 1);

rateSmawk = zeros(nWeights, 2);
rateDp = zeros(nWeights, 2);
rateBisection = zeros(nWeights, 2);
rateMl = zeros(nWeights, 2);
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	% * Threshold by SMAWK
	[thresholdSmawk, dmtcSmawk] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	% * Threshold by dynamic programming
	[thresholdDp, dmtcDp] = threshold_dp(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	% * Threshold by bisection
	[thresholdBisection, dmtcBisection] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	% * Threshold by maximum likelihood detection
	[thresholdMl, dmtcMl] = threshold_ml(equivalentDistribution, receivedPower, symbolRatio);
	% * Compute achievable rate of primary and secondary links
	[~, rateSmawk(iWeight, 1), rateSmawk(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtcSmawk);
	[~, rateDp(iWeight, 1), rateDp(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtcDp);
	[~, rateBisection(iWeight, 1), rateBisection(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtcBisection);
	[~, rateMl(iWeight, 1), rateMl(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtcMl);
end

figure('name', 'Achievable rate region by different threshold design', 'position', [0, 0, 500, 400]);
rateSmawk(end + 1, 1) = max(rateSmawk(:, 1));
rateDp(end + 1, 1) = max(rateDp(:, 1));
rateBisection(end + 1, 1) = max(rateBisection(:, 1));
rateMl(end + 1, 1) = max(rateMl(:, 1));
rateSmawk(end + 2, 2) = max(rateSmawk(:, 2));
rateDp(end + 2, 2) = max(rateDp(:, 2));
rateBisection(end + 2, 2) = max(rateBisection(:, 2));
rateMl(end + 2, 2) = max(rateMl(:, 2));
plotHandle = gobjects(1, 4);
hold all;
plotHandle(1) = plot(rateSmawk(convhull(rateSmawk), 1), rateSmawk(convhull(rateSmawk), 2));
plotHandle(2) = plot(rateDp(convhull(rateDp), 1), rateDp(convhull(rateDp), 2));
plotHandle(3) = plot(rateBisection(convhull(rateBisection), 1), rateBisection(convhull(rateBisection), 2));
plotHandle(4) = plot(rateMl(convhull(rateMl), 1), rateMl(convhull(rateMl), 2));
hold off;
grid minor;
legend('SMAWK', 'Dynamic Programming', 'Bisection', 'Maximum Likelihood', 'Location', 'southwest');
xlabel('Primary rate [nats/s/Hz]');
ylabel('Baskcatter sum rate [nats/channel use]');
xlim([0 inf]);
ylim([0 inf]);
box on;
plot_style(plotHandle);
