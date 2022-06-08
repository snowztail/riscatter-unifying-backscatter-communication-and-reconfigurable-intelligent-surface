clear; setup; cvx_begin; cvx_end; clc; close all; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Result(nVariables, 1) = struct('rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Evaluate rate region for different number of transmit antennas
for iVariable = 1 : nVariables
	% * Set symbol duration ratio
	symbolRatio = Variable.symbolRatio(iVariable);

	% * Generate channels
	directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
	cascadedChannel = zeros(nTxs, nTags);
	for iTag = 1 : nTags
		cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
	end
	equivalentChannel = directChannel + sqrt(scatterRatio) * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

	% * Evaluate rate region
	rate = zeros(2, nWeights);
	distribution = zeros(nStates, nTags, nWeights);
	threshold = zeros(nWeights, nStates ^ nTags + 1);
	beamforming = zeros(nTxs, nWeights);
	for iWeight = 1 : nWeights
		[rate(:, iWeight), distribution(:, :, iWeight), threshold(iWeight, :), beamforming(:, iWeight)] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, weightSet(iWeight), equivalentChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
	end
	Result(iVariable) = struct('rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
end

directory = strcat('data/', mfilename, '/');
data_save;
