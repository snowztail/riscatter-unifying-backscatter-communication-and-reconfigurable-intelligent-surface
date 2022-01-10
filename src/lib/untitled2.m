clear; clc;
nTxs = 1;
nTags = 2;
nStates = 2;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = eps;
reflectRatio = 0.5;
symbolRatio = 10;
noisePower = 1;

constellation = phase_shift_keying(nStates);

precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Primary AWGN channel
indexCombination = nested_combvec(1 : nStates, nTags);
inputCombination = transpose(constellation(indexCombination));
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(reflectRatio) * inputCombination(iInput, :) * cascadedChannel;
end
signalPower = abs(equivalentChannel * precoder) .^ 2;
snr = signalPower / noisePower;

% * Secondary desorted DMTC
[receivedPower, mapIndex] = sort(signalPower + noisePower);
threshold = [0; symbolRatio * movmean(receivedPower, 2, 'Endpoints', 'discard'); Inf];
threshold(2) = threshold(3);
threshold(4) = threshold(3);
threshold_ = threshold;
dmtc = zeros(nInputs, nOutputs);
for iInput = 1 : nInputs
	conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
	for iOutput = 1 : nOutputs
		dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
	end
end
dmtc(mapIndex, mapIndex) = dmtc;

% * Update tag input distribution
[wsr_, primaryRate_, secondaryRate_, inputDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags);
% snr_ = abs(directChannel * precoder) .^ 2 / noisePower;
% primaryRate_ = symbolRatio * log2(1 + snr_);

combinationDistribution = combination_distribution(inputDistribution);
jointDistribution = prod(combinationDistribution, 1);
primaryInformationFunction = primary_information_function(symbolRatio, snr);
primaryRate = jointDistribution * primaryInformationFunction;

% * For a given input distribution, search different threshold
thresholdSet = symbolRatio * nmultichoosek(min(receivedPower) : 5e-2 : max(receivedPower), nStates ^ nTags - 1);
nSets = size(thresholdSet, 1);
secondaryRate = zeros(nSets, 1);
mutualInformation = zeros(nSets, 1);
for iSet = 1 : nSets
	threshold = [0, thresholdSet(iSet, :), Inf];
	dmtc = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
	dmtc(mapIndex, mapIndex) = dmtc;

	secondaryInformationFunction = secondary_information_function(jointDistribution, dmtc);
	secondaryRate(iSet) = jointDistribution * secondaryInformationFunction;
	informationFunction = weight * primaryInformationFunction + (1 - weight) * secondaryInformationFunction;
	mutualInformation(iSet) = jointDistribution * informationFunction;
end
[val, idx] =  max(secondaryRate);
secondaryRate_
val
optThreshold = [0, thresholdSet(idx, :), Inf]

flag = 1;




%% * Local Functions
function [combinationDistribution] = combination_distribution(inputDistribution)
	[nTags, nStates] = size(inputDistribution);
	nInputs = nStates ^ nTags;
	tagSet = transpose(1 : nTags);
	indexCombination = nested_combvec(1 : nStates, nTags);
	combinationDistribution = zeros(nTags, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), tagSet, indexCombination(:, iInput)));
	end
end

function [primaryInformationFunction] = primary_information_function(symbolRatio, snr)
	nInputs = size(snr, 1);
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log2(1 + snr(iInput));
	end
end

function [secondaryInformationFunction] = secondary_information_function(jointDistribution, dmc)
	[nInputs, nOutputs] = size(dmc);
	secondaryInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			secondaryInformationFunction(iInput, iOutput) = dmc(iInput, iOutput) * log2(dmc(iInput, iOutput) / (jointDistribution * dmc(:, iOutput)));
		end
	end
	secondaryInformationFunction = sum(secondaryInformationFunction, 2);
end
