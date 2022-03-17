function [combinationDistribution] = combination_distribution(inputDistribution)
	% Function:
    %	- obtain all combinations of tag input probability
    %
    % Input:
	%	- inputDistribution [nTags * nStates]: input probability distribution
    %
    % Output:
	%	- combinationDistribution [nTags * (nStates ^ nTags)]: each column represents a possible tag input probability combination
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	[nTags, nStates] = size(inputDistribution);
	nInputs = nStates ^ nTags;
	tagSet = transpose(1 : nTags);
	indexCombination = index_combination(nTags, nStates);
	combinationDistribution = zeros(nTags, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), tagSet, indexCombination(:, iInput)));
	end
end
