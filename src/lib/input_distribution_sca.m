function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_sca(nTags, dmtc, weight, symbolRatio, snr, tolerance)
	% Function:
	%	- optimize the tag input distribution by successive convex approximation
    %
    % Input:
	%	- nTags: number of tags
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
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
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
		tolerance = 1e-8;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nInputs = size(dmtc, 1);
	nStates = nthroot(nInputs, nTags);

	% * Initialization
	indexCombination = index_combination(nTags, nStates);
	inputDistribution = ones(nTags, nStates) / nStates;
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc);

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
			weightedSumRateSca = rate_weighted_sum(weight, symbolRatio, snr, auxiliary, dmtc);

			maximize weightedSumRateSca
			subject to
				for iTag = 1 : nTags
					transpose(inputDistribution(iTag, :)) == simplex(nStates);
				end
		cvx_end

		% * Compute actual weighted sum rate
		equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
		weightedSumRate = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc);

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
end
