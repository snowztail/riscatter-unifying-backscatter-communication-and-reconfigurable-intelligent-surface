clear; setup; cvx_begin; cvx_end; clc; config;

% * Initialize struct
Result(nWeights) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

% * Generate channels
directChannel = rxGain * path_loss(frequency, directDistance, directExponent) * fading_ricean(nTxs, 1, directFactor);
cascadedChannel = zeros(nTxs, nTags);
for iTag = 1 : nTags
	cascadedChannel(:, iTag) = rxGain * path_loss(frequency, forwardDistance(iTag), forwardExponent) * fading_ricean(nTxs, 1, forwardFactor) * path_loss(frequency, backwardDistance(iTag), backwardExponent) * fading_ricean(1, 1, backwardFactor);
end
equivalentChannel = directChannel + scatterRatio * cascadedChannel * transpose(constellation(tuple_tag(repmat(transpose(1 : nStates), [1, nTags]))));

% * Clear persistent variables
clear block_coordinate_descent distribution_kkt distribution_cooperation beamforming_pgd;

% * Evaluate rate region
for iWeight = 1 : nWeights
	weight = weightSet(iWeight);
	[rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, cascadedChannel, 'Distribution', 'kkt', 'Beamforming', 'pgd', 'Threshold', 'smawk');
	Result(iWeight) = struct('weight', weight, 'rate', rate, 'distribution', distribution, 'threshold', threshold, 'beamforming', beamforming);
end

%% * Retrieve rate region
rate = zeros(2, nWeights + 3);
for iWeight = 1 : nWeights
	rate(:, iWeight) = Result(iWeight).rate;
end
[rate(1, nWeights + 1), rate(2, nWeights + 2)] = deal(max(rate(1, :)), max(rate(2, :)));
region = rate(:, convhull(transpose(rate)));

%% * Plot rate region
figure('Name', 'Example Primary-(Sum-)Backscatter Rate Region', 'Position', [0, 0, 500, 400]);
object = plot(region(1, :) / log(2), region(2, :) / log(2), 'DisplayName', strcat('Example'));
hold off; legend('Location', 'sw'); grid on; box on; axis tight;
xlabel('Primary Rate [bits/s/Hz]');
ylabel('Total Backscatte Rate [bits/BSP]');
xlim([5, inf]);
style_plot(object);
