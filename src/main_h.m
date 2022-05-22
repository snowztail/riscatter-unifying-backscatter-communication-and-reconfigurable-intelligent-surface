setup; clear; cvx_clear; clc; close all; config;

% * Generate channels
directChannel = sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs));
cascadedChannel = zeros(nTags, nTxs);
for iTag = 1 : nTags
	cascadedChannel(iTag, :) = (sqrt(0.5) * (randn(1, nTxs) + 1i * randn(1, nTxs))) * (sqrt(0.5) * (randn + 1i * randn));
end

% * Obtain equivalent primary channel for each tag state tuple
% [indexTuple, inputTuple] = state_tuple(nTags, nStates);
% indexCombination = index_combination(nTags, nStates);
stateTuple = state_tuple(nTags, nStates);
inputCombination = constellation(stateTuple);
equivalentChannel = zeros(nInputs, nTxs);
for iInput = 1 : nInputs
	equivalentChannel(iInput, :) = directChannel + sqrt(scatterRatio) * inputCombination(iInput, :) * cascadedChannel;
end

% * Evaluate rate regions
rate = zeros(2, nWeights);
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	[rate, inputDistribution, threshold, beamformer] = block_coordinate_descent(weight, nTags, symbolRatio, equivalentChannel, transmitPower, noisePower, nBins, tolerance, 'Input', 'kkt', 'Threshold', 'smawk');
end
