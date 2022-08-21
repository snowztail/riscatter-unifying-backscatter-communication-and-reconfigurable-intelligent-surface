clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = sqrt(path_loss(directDistance, directExponent)) * fading_ricean(nTxs, nRxs, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = sqrt(path_loss(forwardDistance(iTag), forwardExponent)) * fading_ricean(nTxs, nSxs, forwardFactor) * sqrt(path_loss(backwardDistance(iTag), backwardExponent)) * fading_ricean(nSxs, nRxs, backwardFactor);
end

% * Evaluate rate region vs number of reflection states
for iVariable = 1 : nVariables
	% * Set constellation and update equivalent channel
	nStates = Variable(iVariable).nStates;
	constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);
	equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

	% * Clear persistent variables
	clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd threshold_bisection;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, cascadedChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
		Result(iVariable, iWeight) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
