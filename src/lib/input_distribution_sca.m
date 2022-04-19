function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_sca(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc, tolerance)
	% Function:
	%	- optimize the tag input distribution by successive convex approximation
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- nTags: number of tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
	%	- sequentially approximate the equivalent (i.e. product of individual) distribution by its first-order Taylor expansion at previous point
	%	- converges to the KKT points
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 15

	% * Declare default tolerance
	arguments
		weight;
		nTags;
		symbolRatio;
		equivalentChannel;
		noisePower;
		beamformer;
		dmtc;
		tolerance = 1e-6;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nInputs = size(dmtc, 1);
	nStates = nthroot(nInputs, nTags);

	% * Initialize input distribution as uniform distribution
	indexCombination = index_combination(nTags, nStates);
	inputDistribution = ones(nTags, nStates) / nStates;
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);

	% * Iteratively update input distribution by SCA
	isConverged = false;
	while ~isConverged
		inputDistribution_ = inputDistribution;
		weightedSumRate_ = weightedSumRate;
		cvx_begin
			variables inputDistribution(nTags, nStates);
			expressions auxiliary(1, nInputs);
			for iInput = 1 : nInputs
				for iTag = 1 : nTags
					iComplement = setdiff(1 : nTags, iTag);
					auxiliary(iInput) = auxiliary(iInput) + inputDistribution(iTag, indexCombination(iTag, iInput)) * prod(inputDistribution_(sub2ind(size(inputDistribution_), transpose(iComplement), indexCombination(iComplement, iInput))), 1);
				end
				auxiliary(iInput) = auxiliary(iInput) - (nTags - 1) * prod(inputDistribution_(sub2ind(size(inputDistribution_), transpose(1 : nTags), indexCombination(:, iInput))), 1);
			end
			weightedSumRateSca = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, auxiliary, beamformer, dmtc);

			maximize weightedSumRateSca
			subject to
				for iTag = 1 : nTags
					transpose(inputDistribution(iTag, :)) == simplex(nStates);
				end
		cvx_end

		% * Compute actual weighted sum rate
		equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
		weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
end
