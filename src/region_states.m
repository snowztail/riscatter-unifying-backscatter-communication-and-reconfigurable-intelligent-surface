clear; setup; cvx_begin; cvx_end; clc; close all; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, nWeights) = struct('rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Evaluate rate region for different number of reflection states
for iVariable = 1 : nVariables
	% * Set constellation
	nStates = Variable(iVariable).nStates;
	constellation = normalize(qammod(transpose(0 : nStates - 1), nStates), 'norm', Inf);

	% * Generate channels
	directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
	cascadedChannel = zeros(nTxs, nTags);
	for iTag = 1 : nTags
		cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
	end
	equivalentChannel = directChannel + sqrt(scatterRatio) * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

	% * Clear persistent variable
	clear block_coordinate_descent beamforming_pgd;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, weight, equivalentChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
		Result(iVariable, iWeight) = struct('rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
