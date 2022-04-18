function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_exhaustion(weight, nTags, symbolRatio, equivalentChannel, noisePower, precoder, dmtc, resolution)
	% Function:
	%	- obtain the WSR-optimal tag input distribution by (nTags-dimensional) exhaustive search
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- nTags: number of tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- precoder [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- resolution: gap between exhaustive search levels
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
    %	- naive exhaustive search with exponential time complexity
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Declare default resolution
	arguments
		weight;
		nTags;
		symbolRatio;
		equivalentChannel;
		noisePower;
		precoder;
		dmtc;
		resolution = 5e-3;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Initialization
	probabilitySimplex = simplex_probability(nStates, resolution);
	indexCombination = index_combination(nTags, size(probabilitySimplex, 1));
	nCombinations = size(probabilitySimplex, 1) ^ nTags;

	% * Exhaustive search
	inputDistributionSet = cell(nCombinations, 1);
	equivalentDistributionSet = cell(nCombinations, 1);
	weightedSumRateSet = zeros(nCombinations, 1);
	for iCombination = 1 : nCombinations
		inputDistributionSet{iCombination} = probabilitySimplex(indexCombination(:, iCombination), :);
		equivalentDistributionSet{iCombination} = prod(combination_distribution(inputDistributionSet{iCombination}), 1);
		weightedSumRateSet(iCombination) = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistributionSet{iCombination}, precoder, dmtc);
	end
	[weightedSumRate, optimalIndex] = max(weightedSumRateSet);
	inputDistribution = inputDistributionSet{optimalIndex};
	equivalentDistribution = equivalentDistributionSet{optimalIndex};
end
