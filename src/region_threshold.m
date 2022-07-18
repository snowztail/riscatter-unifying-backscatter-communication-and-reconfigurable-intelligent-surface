clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, nWeights) = struct('rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Evaluate rate region by different decision threshold design
for iVariable = 1 : nVariables
	% * Clear persistent variables
	clear block_coordinate_descent;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, cascadedChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', Variable(iVariable).Threshold);
		if any(isnan(rate))
			error('PGD failed. Please ensure DMAC entries are not approaching zero.');
		else
			Result(iVariable, iWeight) = struct('rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
		end
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
