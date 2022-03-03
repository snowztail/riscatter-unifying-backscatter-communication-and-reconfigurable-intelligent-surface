function [rateBound] = rate_bound_randomization(dmtc, subSet, inputDistribution)
	% Function:
	%	- obtain the upper bounds of sum rate on subset of tags
    %
    % Input:
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- subset: the subset of tag indexes
	%	- inputDistribution [nTags * nStates]: input probability distribution
    %
    % Output:
	%	- rateBound: the upper bounds of sum rate on subset of tags
    %
    % Comment:
    %   - Numerically compute the upper bound of the sum rate on the subset for a given input distribution
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 03

	% * Get data
	[nTags, nStates] = size(inputDistribution);
	nOutputs = size(dmtc, 2);

	% * Initialization
	setIndexCombination = combvec_nested(1 : nStates, nTags);
	nSets = size(setIndexCombination, 2);
	setEquivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	subSetEquivalentDistribution = prod(combination_distribution(inputDistribution(subSet, :)), 1);

	% * Obtain upper bounds
	rateBound = zeros(nOutputs, nSets);
	for iOutput = 1 : nOutputs
		for iSet = 1 : nSets
			setIndex = setIndexCombination(:, iSet);
			rateBound(iOutput, iSet) = setEquivalentDistribution(iSet) * dmtc(iSet, iOutput) * log(dmtc(iSet, iOutput) / (subSetEquivalentDistribution * dmtc(all(setIndexCombination(setdiff(1 : nTags, subSet), :) == setIndex(setdiff(1 : nTags, subSet)), 1), iOutput)));
		end
	end
	rateBound = sum(rateBound, 'all');
end
