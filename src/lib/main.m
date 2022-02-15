clear; clc;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = eps;
reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;
nLevels = 2 ^ 6;
constellation = phase_shift_keying(nStates);
confidenceScore = 10;
tolerance = 1e-3;

precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Compute expected received power per primary symbol
indexCombination = nested_combvec(1 : nStates, nTags);
inputCombination = transpose(constellation(indexCombination));
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(reflectRatio) * inputCombination(iInput, :) * cascadedChannel;
end
signalPower = abs(equivalentChannel * precoder) .^ 2;
snr = signalPower / noisePower;
[receivedPower, mapIndex] = sort(signalPower + noisePower);

% * Obtain threshold candidates within empirical interval based on Chebyshev's inequality
lowerBound = symbolRatio * min(receivedPower) - confidenceScore * sqrt(symbolRatio * min(receivedPower));
upperBound = symbolRatio * max(receivedPower) + confidenceScore * sqrt(symbolRatio * max(receivedPower));
if lowerBound <= 0
	thresholdCandidate = [linspace(0, upperBound, nLevels), inf];
else
	thresholdCandidate = [0, linspace(lowerBound, upperBound, nLevels - 1), inf];
end

% * Discretize continuous output and remaps to DMC
dmc = discretize_channel(nInputs, nLevels, symbolRatio, receivedPower, thresholdCandidate);
dmc(mapIndex, :) = dmc;


% * Initialize detection thresholds, DMTC, and input distribution
threshold = [0, thresholdCandidate(discretize(symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'), thresholdCandidate)), inf];
dmtc = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold);
dmtc(mapIndex, mapIndex) = dmtc;
[wsr_, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags);

isConverged = false;
while ~isConverged
	% * Update thresholding scheme
	tic
	[threshold] = thresholding(dmc, equivalentDistribution, thresholdCandidate);
	toc
	tic
	[threshold_smawk] = thresholding_smawk(dmc, equivalentDistribution, thresholdCandidate);
	toc

	% * Update DMTC
	dmtc = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold);
	dmtc(mapIndex, mapIndex) = dmtc;
	dmtc_smawk = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold_smawk);
	dmtc_smawk(mapIndex, mapIndex) = dmtc_smawk;

	% * Update tag input distribution
	[wsr, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags);
	[wsr_smawk, primaryRate_smawk, secondaryRate_smawk, inputDistribution_smawk, equivalentDistribution_smawk] = input_distribution(weight, symbolRatio, snr, dmtc_smawk, nTags);

	% * Test convergence
	isConverged = abs(wsr - wsr_) <= tolerance;
	wsr_ = wsr;
end
