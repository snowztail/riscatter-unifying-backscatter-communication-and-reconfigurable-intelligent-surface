clear; clc;
nTxs = 1;
nTags = 5;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = [eps; 1 - eps];
reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;
nBins = 2 ^ 8;
constellation = psk(nStates);
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
inputDistribution = normr(sqrt(rand(nTags, nStates))) .^ 2;
combinationDistribution = combination_distribution(inputDistribution);
equivalentDistribution = prod(combinationDistribution, 1);
[threshold, dmtc] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
[weightedSumRate, primaryRate, backscatterRate] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc);
% dmtc = normr(sqrt(rand(nStates ^ nTags, nStates ^ nTags))) .^ 2;

% * Update input distribution and detection threshold alternatively
isConverged = false;
while ~isConverged
	[inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_optimization(nTags, dmtc, weight, symbolRatio, snr);
	[inputDistribution1, equivalentDistribution1, weightedSumRate1] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr);
	[threshold, dmtc, backscatterRate] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	[threshold1, dmtc1, backscatterRate1] = threshold_dp(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
	[threshold2, dmtc2, backscatterRate2] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);

	[weightedSumRate_, primaryRate, backscatterRate] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc);
	[weightedSumRate1_, primaryRate1, backscatterRate1] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc1);
	[weightedSumRate2_, primaryRate2, backscatterRate2] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc2);

	isConverged = abs(weightedSumRate_ - weightedSumRate) <= tolerance;
	weightedSumRate = weightedSumRate_;
end