clear; setup; cvx_begin; cvx_end; clc; close all; config;

% * Initialize struct
Result(nWeights) = struct('rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Clear persistent variable
clear block_coordinate_descent beamforming_pgd;

% * Evaluate rate region
for iWeight = 1 : nWeights
	[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weightSet(iWeight), equivalentChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
	Result(iWeight) = struct('rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
end
