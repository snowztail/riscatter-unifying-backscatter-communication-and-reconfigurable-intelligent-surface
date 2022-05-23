setup; clear; cvx_clear; cvx_begin; cvx_end; clc; close all; config;

% * Generate channels
directChannel = path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor, directAoa);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor, forwardAoa(iTag)) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor, backwardAoa(iTag));
end
equivalentChannel = directChannel + sqrt(scatterRatio) * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Evaluate rate region
rate = zeros(2, nWeights);
distribution = zeros(nStates, nTags, nWeights);
threshold = zeros(nWeights, nOutputs + 1);
beamforming = zeros(nTxs, nWeights);
for iWeight = 1 : nWeights
	[rate(:, iWeight), distribution(:, :, iWeight), threshold(iWeight, :), beamforming(:, iWeight)] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, weightSet(iWeight), equivalentChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
	% TODO fine tune tolerance
	if rate(2, iWeight) <= tolerance
		break;
	end
end
