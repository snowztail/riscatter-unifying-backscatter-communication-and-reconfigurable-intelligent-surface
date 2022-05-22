function [beamformer, dmtc, weightedSumRate] = beamformer_sca(weight, symbolRatio, equivalentChannel, transmitPower, noisePower, equivalentDistribution, threshold, tolerance)
	% Function:
	%	- optimize the transmit beamformer to maximize the weighted sum primary-backscatter rate by successive convex approximation
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent primary channel under each tag input combination
	%	- transmitPower: average transmit power at the AP
	%	- noisePower: average noise power at the user
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- threshold [1 * (nOutputs + 1)]: boundaries of decision regions
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
    %	- iteratively update the beamforming matrix until convergence
	%	- extract the best beamformer by Gaussian randomization
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 19

	% * Declare default tolerance
	arguments
		weight;
		symbolRatio;
		equivalentChannel;
		transmitPower;
		noisePower;
		equivalentDistribution;
		threshold;
		tolerance = 1e-6;
	end

	% % ! Use finite threshold to evaluate incomplete gamma function by series representation
	% threshold(end) = 10 * threshold(end - 1);

	% * Get data
	[nInputs, nOutputs] = deal(size(equivalentChannel, 1));
	nTxs = size(equivalentChannel, 2);

	% * Initialize beamformer matrix and DMTC
	beamformer = sqrt(transmitPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
	beamformerMatrix = beamformer * beamformer';
	[dmtc, receivedPower] = channel_discretization(symbolRatio, equivalentChannel, noisePower, beamformer, threshold);
	weightedSumRate = weighted_sum_rate_local(weight, symbolRatio, noisePower, equivalentDistribution, receivedPower, dmtc);

	% * Iteratively update beamforming matrix by SCA
	isConverged = false;
	while ~isConverged
		% * Update iteration index
		weightedSumRate_ = weightedSumRate;
		beamformerMatrix_ = beamformerMatrix;

		% * Solve subproblem
		cvx_begin
		cvx_solver mosek
			variable beamformerMatrix(nTxs, nTxs);
			expression backscatterRateSca;


			for iOutput = 1 : nOutputs
				for iInput = 1 : nInputs
					% * Marginal entropy
					backscatterRateSca = backscatterRateSca + equivalentDistribution(iInput) * (regularizedGamma_(iInput, iOutput) * log(regularizedGamma_(iInput, iOutput)) + (regularizedGamma(iInput, iOutput) - regularizedGamma_(iInput, iOutput)) * (log(regularizedGamma_(iInput, iOutput)) + 1));
				end
				% * Conditional entropy
				backscatterRateSca = backscatterRateSca + entr(equivalentDistribution * regularizedGamma(:, iOutput));
			end
			receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);

			maximize backscatterRateSca
			subject to
				beamformerMatrix == hermitian_semidefinite(nTxs);
				trace(beamformerMatrix) <= transmitPower;
		cvx_end

	end




















































	% * Compute expected received power under all backscatter input combinations and update DMTC
	receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);
	dmtc = dmtc_local(symbolRatio, receivedPower, threshold);

	% * Initialize auxiliary variables associated with regularized incomplete gamma function
	epsilon = epsilon_auxiliary(receivedPower, threshold);
	xi = xi_auxiliary(symbolRatio, epsilon);
	[regularizedGamma, regularizedUpperGamma] = regularized_gamma(epsilon, xi);

	% * Iteratively update beamforming matrix by SCA
	isConverged = false;
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
	while ~isConverged
		% * Update iteration index
		weightedSumRate_ = weightedSumRate;
		[receivedPower_, epsilon_, xi_, regularizedGamma_, regularizedUpperGamma_] = deal(receivedPower, epsilon, xi, regularizedGamma, regularizedUpperGamma);
		[xiLower_, xiUpper_] = deal(xi_);
		[epsilonLower_, epsilonUpper_] = deal(epsilon_);

		cvx_begin
		cvx_solver mosek
			variables beamformerMatrix(nTxs, nTxs) regularizedGamma(nInputs, nOutputs);
			expressions backscatterRateSca;

			for iOutput = 1 : nOutputs
				for iInput = 1 : nInputs
					% * Marginal entropy
					backscatterRateSca = backscatterRateSca + equivalentDistribution(iInput) * (regularizedGamma_(iInput, iOutput) * log(regularizedGamma_(iInput, iOutput)) + (regularizedGamma(iInput, iOutput) - regularizedGamma_(iInput, iOutput)) * (log(regularizedGamma_(iInput, iOutput)) + 1));
				end
				% * Conditional entropy
				backscatterRateSca = backscatterRateSca + entr(equivalentDistribution * regularizedGamma(:, iOutput));
			end
			receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);

			maximize backscatterRateSca
			subject to
				beamformerMatrix == hermitian_semidefinite(nTxs);
				trace(beamformerMatrix) <= transmitPower;
		cvx_end

		% * Update DMTC and compute actual weighted sum rate
		dmtc = dmtc_local(symbolRatio, receivedPower, threshold);
		[weightedSumRate, rate] = rate_weighted_sum_local(weight, symbolRatio, noisePower, equivalentDistribution, receivedPower, dmtc);

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end


end

function [weightedSumRate] = weighted_sum_rate_local(weight, symbolRatio, noisePower, equivalentDistribution, receivedPower, dmtc)
	[nInputs, nOutputs] = deal(size(equivalentDistribution, 2));
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
