function [beamformer, weightedSumRate] = beamformer_bpgd(weight, symbolRatio, equivalentChannel, txPower, noisePower, equivalentDistribution, threshold, tolerance, btls)


	% * Declare default tolerance
	arguments
		weight;
		symbolRatio;
		equivalentChannel;
		txPower;
		noisePower;
		equivalentDistribution;
		threshold;
		% step = 1e-3;
		tolerance = 1e-6;
		btls.Alpha = 0.15;
		btls.Beta = 0.5;
	end

	% ! Use finite threshold to evaluate incomplete gamma function by series representation
	threshold(end) = 10 * threshold(end - 1);

	% * Get data
	% [nInputs, nOutputs] = deal(size(equivalentChannel, 1));
	% nTxs = size(equivalentChannel, 2);

	% * Initialize beamformer
	beamformer = sqrt(txPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
	weightedSumRate = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, threshold);

	% * Projected gradient descent
	isConverged = false;
	while ~isConverged
		% * Update iteration index
		weightedSumRate_ = weightedSumRate;

		% * Compute local gradient
		gradient = gradient_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, threshold, beamformer);

		% ? Backtracking line search
		% step = backtracking_line_search(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, threshold, alpha, beta);
		y1 = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, threshold);
		t = 1;
		beamformerCandidate = beamformer + t * gradient;
		beamformerCandidate = sqrt(txPower) * beamformerCandidate / max(sqrt(txPower), norm(beamformerCandidate));
		y2 = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformerCandidate, threshold);
		% while y2 > y1 + btls.Alpha * t * norm(gradient) ^ 2;
		while y2 < y1 + btls.Alpha * t * norm(gradient) ^ 2
			t = btls.Beta * t;
			beamformerCandidate = beamformer + t * gradient;
			beamformerCandidate = sqrt(txPower) * beamformerCandidate / max(sqrt(txPower), norm(beamformerCandidate));
			y2 = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformerCandidate, threshold);
		end
		beamformer = beamformerCandidate;
		% step = t;
% 		step = 1e-2;

% 		% * Perform unregularized gradient descent
% 		beamformer = beamformer + step * gradient;

% 		% * Project onto subspace (l2-norm ball)
% % 		beamformer = sqrt(txPower) * beamformer / max(1, norm(beamformer));
% 		beamformer = sqrt(txPower) * beamformer / max(sqrt(txPower), norm(beamformer));

		% * Update weighted sum rate
		weightedSumRate = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, threshold);

		% * Test convergence (gradient can be non-zero due to norm constraint)
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance || weightedSumRate < weightedSumRate_;
% 		norm(gradient)
		% weightedSumRate
		% isConverged = norm(gradient) <= 1e-3;
	end
end


function [gradient] = gradient_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, threshold, beamformer)
	[nInputs, nOutputs] = deal(size(equivalentChannel, 1));

	% * sigma and H
	sigma2 = zeros(nInputs, 1);
	% H = zeros(nTxs, nTxs, nInputs);
	H = cell(nInputs, 1);
	for iInput = 1 : nInputs
		% H(:, :, iInput) = equivalentChannel(iInput, :)' * equivalentChannel(iInput, :);
		H{iInput} = equivalentChannel(iInput, :)' * equivalentChannel(iInput, :);
		sigma2(iInput) = beamformer' * H{iInput} * beamformer + noisePower;
	end

	% * u, g and Q
	u = zeros(nInputs, nOutputs + 1);
	g = zeros(nInputs, nOutputs + 1);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs + 1
			u(iInput, iOutput) = threshold(iOutput) / sigma2(iInput);
			g(iInput, iOutput) = exp(-u(iInput, iOutput)) * sum(u(iInput, iOutput) .^ (0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
		end
	end
	Q = - diff(g, 1, 2);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			q(iInput, iOutput) = g(iInput, iOutput) - g(iInput, iOutput + 1);
		end
	end

	% * dg/dw and dQ/dw
	dgdw = cell(nInputs, nOutputs + 1);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs + 1
			% dQdw(iInput, iOutput) = H{iInput} * beamformer / sigma2(iInput) ^ 2 * (threshold(iOutput + 1) * exp(-u(iInput, iOutput + 1)) * (sum(((1 : symbolRatio - 1) - u(iInput, iOutput + 1)) * u(iInput, iOutput + 1)) - 1) )
			dgdw{iInput, iOutput} = - threshold(iOutput) * H{iInput} * beamformer / sigma2(iInput) ^ 2 * exp(-u(iInput, iOutput)) * (sum(((1 : symbolRatio - 1) - u(iInput, iOutput)) .* u(iInput, iOutput) .^ (0 : symbolRatio - 2) ./ factorial(1 : symbolRatio - 1)) - 1);
		end
	end
	dQdw = cell(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			dQdw{iInput, iOutput} = dgdw{iInput, iOutput} - dgdw{iInput, iOutput + 1};
		end
	end

	% * dI/dw
	dIdw = 0;
	for iInput = 1 : nInputs
		dIdw = dIdw + weight * equivalentDistribution(iInput) * symbolRatio * H{iInput} * beamformer / sigma2(iInput);
		for iOutput = 1 : nOutputs
			dIdw = dIdw + (1 - weight) * equivalentDistribution(iInput) * ((log(Q(iInput, iOutput) / (equivalentDistribution * Q(:, iOutput))) + 1) * dQdw{iInput, iOutput} - Q(iInput, iOutput) * sum(equivalentDistribution .* cat(2, dQdw{:, iOutput}), 2) / (equivalentDistribution * Q(:, iOutput)));
		end
	end

	gradient = 2 * dIdw;
end

function [dmtc] = dmtc_local(symbolRatio, receivedPower, threshold)
	nInputs = size(receivedPower, 1);
	nOutputs = size(threshold, 2) - 1;
	dmtc = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
end

function [weightedSumRate] = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, threshold)
	[nInputs, nOutputs] = deal(size(equivalentDistribution, 2));
	receivedPower = abs(equivalentChannel * beamformer) .^ 2 + noisePower;
	dmtc = dmtc_local(symbolRatio, receivedPower, threshold);
	primaryRate = equivalentDistribution * symbolRatio * log(1 + receivedPower / noisePower);
	backscatterInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	backscatterRate = equivalentDistribution * sum(backscatterInformationFunction, 2);
	weightedSumRate = [weight, 1 - weight] * [primaryRate; backscatterRate];
end
