clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'region_')));

% * Initialize struct
Reference(1, nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);
Result(nVariables, nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = sqrt(path_loss(directDistance, directExponent)) * fading_ricean(nTxs, nRxs, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = sqrt(path_loss(forwardDistance(iTag), forwardExponent)) * fading_ricean(nTxs, nSxs, forwardFactor) * sqrt(path_loss(backwardDistance(iTag), backwardExponent)) * fading_ricean(nSxs, nRxs, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Evaluate rate region vs cascaded channel estimation error
for iVariable = 1 : nVariables
	% * Set cascaded channel estimation error
	channelVariance = Variable(iVariable).channelVariance;
	cascadedChannelEstimation = cascadedChannel .* (ones(size(cascadedChannel)) + sqrt(0.5 * channelVariance) * (randn(size(cascadedChannel)) + 1i * randn(size(cascadedChannel))));
	equivalentChannelEstimation = directChannel + scatterRatio * cascadedChannelEstimation * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

	% * Clear persistent variables
	clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd threshold_bisection;

	% * Evaluate rate region
	for iWeight = 1 : nWeights
		weight = weightSet(iWeight);
		[~, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannelEstimation, cascadedChannelEstimation, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');

		signalPower = abs(equivalentChannel' * beamforming) .^ 2;
		receivePower = signalPower + noisePower;
		snr = signalPower / noisePower;
		dmac = dmc_integration(symbolRatio, receivePower, threshold);
		[~, sortIndex] = sort(receivePower);
		dmac(:, sortIndex) = dmac;
		[~, rate] = rate_weighted(snr, equivalentDistribution, dmac, weight);

		Result(iVariable, iWeight) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
	end
end

directory = strcat('data/', mfilename, '/');
data_save;
