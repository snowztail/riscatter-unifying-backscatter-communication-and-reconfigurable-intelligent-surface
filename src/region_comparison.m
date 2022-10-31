clear; setup; cvx_begin; cvx_end; clc; config_comparison; 
close all;

% * Initialize struct
ResultRiscatter(nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = sqrt(path_loss(directDistance, directExponent)) * fading_ricean(nTxs, nRxs, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = sqrt(path_loss(forwardDistance(iTag), forwardExponent)) * fading_ricean(nTxs, nSxs, forwardFactor) * sqrt(path_loss(backwardDistance(iTag), backwardExponent)) * fading_ricean(nSxs, nRxs, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Clear persistent variables
clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd threshold_bisection;

% * Obtain achievable rates for legacy, BBC, AmBC, SR, RIS and RIScatter
[rateLegacy] = benchmark_legacy(transmitPower, noisePower, directChannel);
[rateBbc, distributionBbc, thresholdBbc] = benchmark_bbc(nTags, symbolRatio, transmitPower, noisePower, equivalentChannel);
[rateAmbc, distributionAmbc, thresholdAmbc] = benchmark_ambc(nTags, symbolRatio, transmitPower, noisePower, equivalentChannel, scatterRatio, directChannel, cascadedChannel);
[rateSr, distributionSr] = benchmark_sr(nTags, transmitPower, noisePower, equivalentChannel);
[rateRis, distributionRis] = benchmark_ris(nTags, transmitPower, noisePower, equivalentChannel);
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	[rateRiscatter, distributionRiscatter, thresholdRiscatter, beamformingRiscatter] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, cascadedChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
	ResultRiscatter(iWeight) = struct('weight', weight, 'rate', rateRiscatter, 'distribution', distributionRiscatter, 'threshold', thresholdRiscatter, 'beamforming', beamformingRiscatter);
end

save(strcat('data/', mfilename));
