clear; cvx_clear; clc; setup;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = eps;
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
[weightedSumRate, primaryRate, backscatterRate] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc);

% * Update input distribution and detection threshold alternatively
isConverged = false;
while ~isConverged
	% * Joint input with full transmit cooperation
	[jointDistribution, equivalentDistributionCooperation, weightedSumRateCooperation] = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
	% * Input distribution by exhaustive search
	[inputDistributionExhaustion, equivalentDistributionExhaustion, weightedSumRateExhaustion] = input_distribution_exhaustion(nTags, dmtc, weight, symbolRatio, snr);
	% * Input distribution by SCA
	% [inputDistributionSca, equivalentDistributionSca, weightedSumRateSca] = input_distribution_sca(nTags, dmtc, weight, symbolRatio, snr);
	% * Input distribution by KKT solution
	[inputDistributionKkt, equivalentDistributionKkt, weightedSumRateKkt] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr);
	% * Individual input recovery by randomization, marginalization, and decomposition
	[inputDistributionRandomization, equivalentDistributionRandomization, weightedSumRateRandomization] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionMarginalization, equivalentDistributionMarginalization, weightedSumRateMarginalization] = recovery_marginalization(jointDistribution, dmtc, weight, symbolRatio, snr);
	[inputDistributionDecomposition, equivalentDistributionDecomposition, weightedSumRateDecomposition] = recovery_decomposition(jointDistribution, dmtc, weight, symbolRatio, snr);
	% * Thresholding
	[thresholdSmawk, dmtcSmawk, backscatterRateSmawk] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	[thresholdDp, dmtcDp, backscatterRateDp] = threshold_dp(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	[thresholdBisection, dmtcBisection, backscatterRateBisection] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	[thresholdMl, dmtcMl, backscatterRateMl] = threshold_ml(equivalentDistribution, receivedPower, symbolRatio);
end
