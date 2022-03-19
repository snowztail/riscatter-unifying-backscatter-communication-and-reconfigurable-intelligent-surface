function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_exhaustion(nTags, dmtc, weight, symbolRatio, snr, resolution)
	% Function:
	%	- obtain the WSR-optimal tag input distribution by (nTags-dimensional) exhaustive search
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
    %	- naive exhaustive search with exponential time complexity
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Declare default resolution
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
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
		weightedSumRateSet(iCombination) = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistributionSet{iCombination}, dmtc);
	end
	[weightedSumRate, optimalIndex] = max(weightedSumRateSet);
	inputDistribution = inputDistributionSet{optimalIndex};
	equivalentDistribution = equivalentDistributionSet{optimalIndex};
end
