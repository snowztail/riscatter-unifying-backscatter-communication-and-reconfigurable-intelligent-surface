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
[threshold, dmtc] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);

rateCooperation = zeros(nWeights, 2);
rateExhaustion = zeros(nWeights, 2);
rateKkt = zeros(nWeights, 2);
rateRandomization = zeros(nWeights, 2);
rateMarginalization = zeros(nWeights, 2);
rateDecomposition = zeros(nWeights, 2);
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	% * Joint input with full transmit cooperation
	[jointDistribution, equivalentDistributionCooperation, weightedSumRateUpperBound] = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
	% * Input distribution by exhaustive search
	[inputDistributionExhaustion, equivalentDistributionExhaustion, weightedSumRateExhaustion] = input_distribution_exhaustion(nTags, dmtc, weight, symbolRatio, snr);
	% * Input distribution by KKT solution
	[inputDistributionKkt, equivalentDistributionKkt, weightedSumRateKkt] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr);
	% * Individual input recovery by randomization, marginalization, and decomposition
	[inputDistributionRandomization, equivalentDistributionRandomization, weightedSumRateRandomization] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionMarginalization, equivalentDistributionMarginalization, weightedSumRateMarginalization] = recovery_marginalization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionDecomposition, equivalentDistributionDecomposition, weightedSumRateDecomposition] = recovery_decomposition(jointDistribution, dmtc, weight, symbolRatio, snr);
	% * Compute achievable rate of primary and secondary links
	[~, rateCooperation(iWeight, 1), rateCooperation(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionCooperation, dmtc);
	[~, rateExhaustion(iWeight, 1), rateExhaustion(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionExhaustion, dmtc);
	[~, rateKkt(iWeight, 1), rateKkt(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionKkt, dmtc);
	[~, rateRandomization(iWeight, 1), rateRandomization(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionRandomization, dmtc);
	[~, rateMarginalization(iWeight, 1), rateMarginalization(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionMarginalization, dmtc);
	[~, rateDecomposition(iWeight, 1), rateDecomposition(iWeight, 2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionDecomposition, dmtc);
end

figure('name', 'Achievable rate region by different input distribution design', 'position', [0, 0, 500, 400]);
rateCooperation(end + 2, 2) = max(rateCooperation(:, 2));
rateExhaustion(end + 2, 2) = max(rateExhaustion(:, 2));
rateKkt(end + 2, 2) = max(rateKkt(:, 2));
rateRandomization(end + 2, 2) = max(rateRandomization(:, 2));
rateMarginalization(end + 2, 2) = max(rateMarginalization(:, 2));
rateDecomposition(end + 2, 2) = max(rateDecomposition(:, 2));
plotHandle = gobjects(1, 6);
hold all;
plotHandle(1) = plot(rateCooperation(convhull(rateCooperation), 1), rateCooperation(convhull(rateCooperation), 2));
plotHandle(2) = plot(rateExhaustion(convhull(rateExhaustion), 1), rateExhaustion(convhull(rateExhaustion), 2));
plotHandle(3) = plot(rateKkt(convhull(rateKkt), 1), rateKkt(convhull(rateKkt), 2));
plotHandle(4) = plot(rateRandomization(convhull(rateRandomization), 1), rateRandomization(convhull(rateRandomization), 2));
plotHandle(5) = plot(rateMarginalization(convhull(rateMarginalization), 1), rateMarginalization(convhull(rateMarginalization), 2));
plotHandle(6) = plot(rateDecomposition(convhull(rateDecomposition), 1), rateDecomposition(convhull(rateDecomposition), 2));
hold off;
grid minor;
legend('Tag Cooperation', 'Exhaustive Search', 'KKT', 'Randomization', 'Marginalization', 'Decomposition', 'Location', 'southwest');
xlabel('Primary rate [nats/s/Hz]');
ylabel('Baskcatter sum rate [nats/channel use]');
xlim([0 inf]);
ylim([0 inf]);
box on;
plot_style(plotHandle);
