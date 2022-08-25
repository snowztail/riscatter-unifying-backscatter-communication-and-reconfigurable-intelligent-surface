clear; setup; cvx_begin; cvx_end; clc; run(strcat('config_', erase(mfilename, 'rate_')));

% * Initialize struct
Result(nVariables, 1) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = sqrt(path_loss(directDistance, directExponent)) * fading_ricean(nTxs, nRxs, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = sqrt(path_loss(forwardDistance(iTag), forwardExponent)) * fading_ricean(nTxs, nSxs, forwardFactor) * sqrt(path_loss(backwardDistance(iTag), backwardExponent)) * fading_ricean(nSxs, nRxs, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Fix input distribution and active beamforming
distribution = normalize(ones(nStates, nTags), 'norm', 1);
beamforming = sqrt(transmitPower) * directChannel / norm(directChannel);
equivalentDistribution = prod(tuple_tag(distribution), 2);
receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
snr = receivePower / noisePower;

% * Evaluate achievable rates vs number of output quantization bits
for iVariable = 1 : nVariables
	% * Set quantization bins
	nBins = 2 ^ Variable(iVariable).nBits;

	% * Update thresholds
	thresholdDomain = domain_threshold(symbolRatio, nBins, receivePower);
	binDmc = dmc_integration(symbolRatio, receivePower, thresholdDomain);
	threshold = threshold_smawk(equivalentDistribution, thresholdDomain, binDmc);

	% * Construct DMTC and recover i/o mapping
	dmac = dmc_integration(symbolRatio, receivePower, threshold);
	[~, sortIndex] = sort(receivePower);
	dmac(:, sortIndex) = dmac;
	[wsr, rate] = rate_weighted(weight, snr, equivalentDistribution, dmac);
	Result(iVariable) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
end

directory = strcat('data/', mfilename, '/');
data_save;
