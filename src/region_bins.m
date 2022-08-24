clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = sqrt(path_loss(directDistance, directExponent)) * fading_ricean(nTxs, nRxs, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = sqrt(path_loss(forwardDistance(iTag), forwardExponent)) * fading_ricean(nTxs, nSxs, forwardFactor) * sqrt(path_loss(backwardDistance(iTag), backwardExponent)) * fading_ricean(nSxs, nRxs, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Evaluate rate region vs number of output quantization bits
for iVariable = 1 : nVariables
	% * Set quantization bins
	nBins = Variable(iVariable).nBins;

	% * Clear persistent variables
	clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd threshold_bisection;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
		Result(iVariable, iWeight) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
