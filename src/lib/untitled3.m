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
tolerance = eps;

precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Typical primary symbol power
indexCombination = nested_combvec(1 : nStates, nTags);
inputCombination = transpose(constellation(indexCombination));
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(reflectRatio) * inputCombination(iInput, :) * cascadedChannel;
end
signalPower = abs(equivalentChannel * precoder) .^ 2;
snr = signalPower / noisePower;
[receivedPower, mapIndex] = sort(signalPower + noisePower);

% * Quantize output power per secondary symbol for threshold design
% Some Erlang tables can be used here to improve accuracy
% quantizedPower = linspace(0, max(receivedPower) + 3 * sqrt(symbolRatio * max(receivedPower) ^ 2), nLevels + 1);
% quantizedPower = [0, linspace(symbolRatio * min(receivedPower), symbolRatio * max(receivedPower), nLevels - 1), inf];
quantizedPower = [0, linspace(0.25 * symbolRatio * min(receivedPower), 2.5 * symbolRatio * max(receivedPower), nLevels - 1), inf];

% * Quantize continuous-output channel to construct and desort DMC
dmc = discretize_channel(nInputs, nLevels, symbolRatio, receivedPower, quantizedPower);
dmc(mapIndex, :) = dmc;


% * Initialize detection thresholds, DMTC, and input distribution
threshold = [0, quantizedPower(discretize(symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'), quantizedPower)), inf];
dmtc = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold);
dmtc(mapIndex, mapIndex) = dmtc;
[wsr_, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags);

isConverged = false;
while ~isConverged
	% * Update thresholding scheme
	tic
	[threshold] = thresholding(dmc, equivalentDistribution, quantizedPower);
	toc
	tic
	[threshold_smawk] = thresholding_smawk(dmc, equivalentDistribution, quantizedPower);
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
