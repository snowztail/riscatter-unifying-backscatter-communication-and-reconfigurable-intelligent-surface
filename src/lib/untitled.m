clear; clc;
nTxs = 1;
nTags = 2;
nStates = 4;
[nInputs, nOutputs] = deal(nStates ^ nTags);
weight = 1 - eps;
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
dmtc = zeros(nInputs, nOutputs);
for iInput = 1 : nInputs
	conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
	for iOutput = 1 : nOutputs
		dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
	end
end
dmtc(mapIndex, mapIndex) = dmtc;

% * Update tag input distribution
% [capacity, inputDistribution] = blahut_arimoto(dmtc);
% [capacity, inputDistribution] = multiuser_blahut_arimoto(dmtc, nTags);
[wsr, primaryRate, secondaryRate, inputDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags);
snr_ = abs(directChannel * precoder) .^ 2 / noisePower;
primaryRate_ = symbolRatio * log2(1 + snr_);
