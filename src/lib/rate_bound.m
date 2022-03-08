function [rateBound] = rate_bound(dmtc, subSet, jointArray, complementArray)
	% Function:
	%	- obtain the upper bounds of sum rate on subset of tags
    %
    % Input:
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- subset: the subset of tag indexes
	%	- jointArray [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
	%	- complementArray [nStates * ... (cardinate of subset complement) ... * nStates]: the joint input probability matrix of tags that belongs to subset complement
    %
    % Output:
	%	- rateBound: the upper bounds of sum rate on subset of tags
    %
    % Comment:
    %	- Reverse the order of dimensions of the input probability matrices for the ease of implementation
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Get data
	nTags = ndims_modified(jointArray);
	nStates = size(jointArray, 1);
	nOutputs = size(dmtc, 2);

	% * Initialization
	complementSet = setdiff(1 : nTags, subSet);
	setIndexCombination = combvec_nested(1 : nStates, nTags);
	subSetIndexCombination = combvec_nested(1 : nStates, length(subSet));
	complementSetIndexCombination = combvec_nested(1 : nStates, length(complementSet));
	nSubSets = size(subSetIndexCombination, 2);
	nComplementSets = size(complementSetIndexCombination, 2);

	% * Obtain upper bounds
	marginalDistribution = cvx(zeros(nOutputs, nComplementSets, nSubSets));
	conditionalEntropy = cvx(zeros(nOutputs, 1));
	relativeEntropy = cvx(zeros(nOutputs, nComplementSets));
	for iOutput = 1 : nOutputs
		% * Simplify operation by rearranging array dimensions in natural order
		conditionalEntropy(iOutput) = transpose(vec(permute(jointArray, ndims_modified(jointArray) : -1 : 1))) * entr(dmtc(:, iOutput));
		for iComplementSet = 1 : nComplementSets
			complementSetIndex = transpose(complementSetIndexCombination(:, iComplementSet));
			complementSetIndexCell = num2cell(complementSetIndex);
			for iSubSet = 1 : nSubSets
				subSetIndex = transpose(subSetIndexCombination(:, iSubSet));
				setIndex = [subSetIndex complementSetIndex];
				setIndex([subSet complementSet]) = setIndex;
				setIndexCell = num2cell(setIndex);
				marginalDistribution(iOutput, iComplementSet, iSubSet) = jointArray(setIndexCell{:}) * dmtc(all(setIndexCombination == transpose(setIndex), 1), iOutput);
			end
			if ~isempty(complementArray(complementSetIndexCell{:}))
				relativeEntropy(iOutput, iComplementSet) = rel_entr(sum(marginalDistribution(iOutput, iComplementSet, :)), complementArray(complementSetIndexCell{:}));
			else
				relativeEntropy(iOutput, iComplementSet) = - entr(sum(marginalDistribution(iOutput, iComplementSet, :)));
			end
		end
	end
	rateBound = - sum(conditionalEntropy) - sum(sum(relativeEntropy));
end
