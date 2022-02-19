function [combinationDistribution] = combination_distribution(inputDistribution)
	% Function:
    %   - obtain all combinations of tag input probability
    %
    % Input:
	%	- nTags: number of tags
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	[nTags, nStates] = size(inputDistribution);
	nInputs = nStates ^ nTags;
	tagSet = transpose(1 : nTags);
	indexCombination = combvec_nested(1 : nStates, nTags);
	combinationDistribution = zeros(nTags, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), tagSet, indexCombination(:, iInput)));
	end
end
