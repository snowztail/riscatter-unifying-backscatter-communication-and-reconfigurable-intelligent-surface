clear; clc;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = [eps; 1 - eps];
reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;
nBins = 2 ^ 8;
constellation = psk(nStates);
confidenceScore = 10;
tolerance = 1e-3;
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

% * Initialize input distribution
inputDistribution = normr(sqrt(rand(nTags, nStates))) .^ 2;
combinationDistribution = combination_distribution(inputDistribution);
equivalentDistribution = prod(combinationDistribution, 1);
[threshold, backscatterMutualInformation] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);

% * Initialize detection thresholds, DMTC, and input distribution
threshold = [0, thresholdCandidate(discretize(symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'), thresholdCandidate)), inf];
% threshold2 = thresholdCandidate(discretize(symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'), thresholdCandidate));
% threshold1 = [0, transpose(interp1(thresholdCandidate(1 : end - 1), thresholdCandidate(1 : end - 1), symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'), 'nearest')), inf];
dmtc_ = discretize_channel(symbolRatio, receivedPower, threshold);
dmtc_(mapIndex, mapIndex) = dmtc_;
[wsr_, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc_, nTags, nStates);
isConverged = false;
while ~isConverged
	% * Eliminate unused letters
	inputIndex = find(equivalentDistribution >= tolerance);
	nOutputs = length(inputIndex);
	p(mapIndex) = receivedPower;
	[receivedPower1, mapIndex1] = sort(signalPower(inputIndex) + noisePower);


	dmtc1 = dmtc_(inputIndex, inputIndex);
	dmtc1 = dmtc1 ./ sum(dmtc1, 2);
% 	[wsr0, primaryRate0, secondaryRate0, inputDistribution0, equivalentDistribution0] = input_distribution(weight, symbolRatio, snr(inputIndex), dmtc1, 1, nStates);
	dmtc1(4, 4) = 0;
	[wsr1, primaryRate1, secondaryRate1, inputDistribution1, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc1, nTags, nStates);

	% * Update thresholding scheme
	[threshold4] = thresholding_smawk(dmc(inputIndex, :), equivalentDistribution(inputIndex), thresholdCandidate);
	dmtc4 = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold4);
	dmtc4(mapIndex, mapIndex1) = dmtc4;
	dmtc4(dmtc4 < eps) = eps;
	secondaryInformationFunction = secondary_information_function(equivalentDistribution, dmtc4);
	mutualInformation4 = equivalentDistribution * secondaryInformationFunction;

% 	[threshold_] = thresholding_smawk_1(dmc, equivalentDistribution, nOutputs, thresholdCandidate);
	[threshold3] = thresholding_smawk(dmc, equivalentDistribution, thresholdCandidate);

	[threshold_smawk, mutualInformation_smawk] = thresholding_smawk_1(dmc, equivalentDistribution, thresholdCandidate, symbolRatio, receivedPower, mapIndex);
	[threshold_bis, mutualInformation_bis] = thresholding_bisection(dmc, equivalentDistribution, thresholdCandidate, symbolRatio, receivedPower, mapIndex);

% 	[threshold_smawk4, mutualInformation_smawk4] = thresholding_smawk_1(dmc(inputIndex, :), equivalentDistribution(inputIndex), thresholdCandidate, symbolRatio, receivedPower, mapIndex);
	% * Update DMTC
% 	dmtc = discretize_channel(nOutputs, nOutputs, symbolRatio, receivedPower1, threshold);
	dmtc = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold);
% 	dmtc(mapIndex1, mapIndex1) = dmtc;
	dmtc(mapIndex, mapIndex1) = dmtc;

	dmtc3 = discretize_channel(nInputs, nInputs, symbolRatio, receivedPower, threshold3);
	dmtc3(mapIndex, mapIndex) = dmtc3;

% 	primaryInformationFunction = information_function_primary(symbolRatio, snr(inputIndex));
% 	secondaryInformationFunction = secondary_information_function(equivalentDistribution(inputIndex), dmtc);
% 	primaryRate = equivalentDistribution(inputIndex) * primaryInformationFunction;
% 	secondaryRate = equivalentDistribution(inputIndex) * secondaryInformationFunction;
% 	informationFunction = weight * primaryInformationFunction + (1 - weight) * secondaryInformationFunction;
% 	mutualInformation = equivalentDistribution(inputIndex) * informationFunction;

	% * Update tag input distribution
	[wsr, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags, nStates);
	dmtc2 = dmtc;
	dmtc2(4,4) = 0;
	[wsr2, primaryRate2, secondaryRate2, inputDistribution2, equivalentDistribution2] = input_distribution(weight, symbolRatio, snr, dmtc2, nTags, nStates);

	[wsr3, primaryRate3, secondaryRate3, inputDistribution3, equivalentDistribution3] = input_distribution(weight, symbolRatio, snr, dmtc3, nTags, nStates);

	% * Test convergence
	isConverged = abs(wsr - wsr_) <= tolerance;
	wsr_ = wsr;
end

function [secondaryInformationFunction] = secondary_information_function(equivalentDistribution, dmtc)
	[nInputs, nOutputs] = size(dmtc);
	secondaryInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			secondaryInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log2(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	secondaryInformationFunction = sum(secondaryInformationFunction, 2);
end
