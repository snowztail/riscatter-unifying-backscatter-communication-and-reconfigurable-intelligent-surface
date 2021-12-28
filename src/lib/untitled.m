clear; clc;
nTxs = 1;
nTags = 2;
nStates = 4;
[nInputs, nOutputs] = deal(nStates ^ nTags);
reflectRatio = 0.5;
nRatio = 10;
constellation = phase_shift_keying(nStates);
noiseVariance = 0;

precoder = normc(randn(nTxs, 1));
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

indexCombination = nested_combvec(1 : nStates, nTags);
inputCombination = transpose(constellation(indexCombination));
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(reflectRatio) * inputCombination(iInput, :) * cascadedChannel;
end

% * Sorted DMTC
[powerLevel, mapIndex] = sort(abs(equivalentChannel * precoder) .^ 2);
threshold = [0; nRatio * movmean(powerLevel, 2, 'Endpoints', 'discard'); Inf];
dmtc = zeros(nInputs, nOutputs);
for iInput = 1 : nInputs
	conditionalEnergy = @(z) (z .^ (nRatio - 1) .* exp(-z ./ powerLevel(iInput))) ./ (powerLevel(iInput) .^ nRatio .* gamma(nRatio));
	for iOutput = 1 : nOutputs
		dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
	end
end

% * Desorted DMTC
dmtc(mapIndex, mapIndex) = dmtc;

% [capacity, inputDistribution] = blahut_arimoto(dmtc);
[capacity, inputDistribution] = multiuser_blahut_arimoto(dmtc, nTags);
