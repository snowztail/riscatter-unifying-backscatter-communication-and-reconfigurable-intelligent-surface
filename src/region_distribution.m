clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Evaluate rate region by different input distribution design
for iVariable = 1 : nVariables
	% * Clear persistent variables
	clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, cascadedChannel, 'Distribution', Variable(iVariable).Distribution, 'Beamforming', 'pgd', 'Threshold', 'smawk', 'Recovery', Variable(iVariable).Recovery);
		Result(iVariable, iWeight) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
