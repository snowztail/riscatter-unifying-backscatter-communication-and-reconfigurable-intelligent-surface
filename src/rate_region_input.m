clear; cvx_clear; clc; close all; setup;
nTxs = 1;
nTags = 3;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
% weight = [eps; 1 - eps];

nWeights = 2e1;
% weightSet = [linspace(0, 1, nWeights); linspace(1, 0, nWeights)];
% weightSet = [linspace(1 - eps, eps, nWeights); linspace(eps, 1 - eps, nWeights)];
% weightSet = [1 - exp(-(linspace(0, 10, nWeights))); exp(-(linspace(0, 8, nWeights)))];
% weightSet = [logspace(log(eps), log(1 - eps), nWeights); 1 - logspace(log(eps), log(1 - eps), nWeights)];
% weightSet = [tanh(0 : nWeights - 1); 1 - tanh(0 : nWeights - 1)];
% weight = [exp(-(linspace(0, pi, nWeights) .^ 2)); 1 - exp(-(linspace(0, pi, nWeights) .^ 2))];
weightSet = [linspace(0, 0.2, nWeights - 1), 1 - eps; linspace(1, 0.8, nWeights - 1), eps];

reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;
nBins = 2 ^ 8;
constellation = qammod(0 : nStates - 1, nStates);
confidenceScore = 10;
tolerance = eps;
precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Compute expected received power per primary symbol
indexCombination = combvec_nested(1 : nStates, nTags);
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
dmc = discretize_channel(thresholdCandidate, receivedPower, symbolRatio);

% * Initialize input distribution, detection threshold, and DMTC
inputDistribution = rand_normalized([nTags, nStates], 2);
combinationDistribution = combination_distribution(inputDistribution);
equivalentDistribution = prod(combinationDistribution, 1);
[threshold, dmtc] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);

rateJoint = zeros(nWeights, 2);
rateRandomization = zeros(nWeights, 2);
rateMarginalization = zeros(nWeights, 2);
rateDecomposition = zeros(nWeights, 2);
rateKkt = zeros(nWeights, 2);
for iWeight = 1 : nWeights
	% * Retrieve weight
	weight = weightSet(:, iWeight);

	% * Joint input optimization
	[jointDistribution, equivalentDistributionJoint, weightedSumRateUpperBound] = input_distribution_optimization(nTags, dmtc, weight, symbolRatio, snr);

	% * Individual input recovery by randomization, marginalization, and decomposition
	[inputDistributionRandomization, equivalentDistributionRandomization, weightedSumRateRandomization] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionMarginalization, equivalentDistributionMarginalization, weightedSumRateMarginalization] = recovery_marginalization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionDecomposition, equivalentDistributionDecomposition, weightedSumRateDecomposition] = recovery_decomposition(jointDistribution, dmtc, weight, symbolRatio, snr);

	% * Input distribution by KKT solution
	[inputDistributionKkt, equivalentDistributionKkt, weightedSumRateKkt] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr);

	% * Compute achievable rate of primary and secondary links
	[~, rateJoint(iWeight, 1), rateJoint(iWeight, 2)] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionJoint, dmtc);
	[~, rateRandomization(iWeight, 1), rateRandomization(iWeight, 2)] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionRandomization, dmtc);
	[~, rateMarginalization(iWeight, 1), rateMarginalization(iWeight, 2)] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionMarginalization, dmtc);
	[~, rateDecomposition(iWeight, 1), rateDecomposition(iWeight, 2)] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionDecomposition, dmtc);
	[~, rateKkt(iWeight, 1), rateKkt(iWeight, 2)] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionKkt, dmtc);
end

figure('name', 'Achievable rate region by different input distribution design', 'position', [0, 0, 500, 400]);
rateJoint(end + 2, 2) = max(rateJoint(:, 2));
rateRandomization(end + 2, 2) = max(rateRandomization(:, 2));
rateMarginalization(end + 2, 2) = max(rateMarginalization(:, 2));
rateDecomposition(end + 2, 2) = max(rateDecomposition(:, 2));
rateKkt(end + 2, 2) = max(rateKkt(:, 2));
hold all;
plot(rateJoint(convhull(rateJoint), 1), rateJoint(convhull(rateJoint), 2));
plot(rateRandomization(convhull(rateRandomization), 1), rateRandomization(convhull(rateRandomization), 2));
plot(rateMarginalization(convhull(rateMarginalization), 1), rateMarginalization(convhull(rateMarginalization), 2));
plot(rateDecomposition(convhull(rateDecomposition), 1), rateDecomposition(convhull(rateDecomposition), 2));
plot(rateKkt(convhull(rateKkt), 1), rateKkt(convhull(rateKkt), 2));
hold off;
grid minor;
legend('Outer Bound', 'Randomization', 'Marginalization', 'Decomposition', 'KKT');
xlabel('Primary rate [nats/s/Hz]');
ylabel('Baskcatter sum rate [nats/channel use]');
xlim([0 inf]);
ylim([0 inf]);
box on;
