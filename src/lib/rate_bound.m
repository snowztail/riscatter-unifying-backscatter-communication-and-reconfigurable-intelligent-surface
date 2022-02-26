function [rateBound] = rate_bound(dmtc, subset, jointArray, complementArray)
	% Function:
	%	- obtain the upper bounds of sum rate on subset of tags
    %
    % Input:
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- jointArray [nStates * ... (nTags) ... * nStates]: the joint input probability matrix of all tags
	%	- complementArray [nStates * ... (cardinate of subset complement) ... * nStates]: the joint input probability matrix of tags that belongs to subset complement
    %
    % Output:
	%	- rateBound: the upper bounds of sum rate on subset of tags
    %
    % Comment:
    %   - Reverse the order of dimensions of the input probability matrices for the ease of implementation
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Get data
	nTags = ndims(jointArray);
	nStates = size(jointArray, 1);
	nOutputs = size(dmtc, 2);

	% * Initialization
	indexCombination = combvec_nested(1 : nStates, nTags);

	% * Permute dimensions to align with natural order (indexCombination)
	P = permute(jointArray, [ndims(jointArray) : -1 : 1]);
	P_c = permute(complementArray, [ndims(complementArray) : -1 : 1]);

	% * Formulate upper bounds
	conditionalEntropy = zeros(nOutputs, 1);
	for iOutput = 1 : nOutputs
		conditionalEntropy(iOutput) = transpose(vec(P)) * entr(dmtc(:, iOutput));
% 			ft1(iOutput) = transpose(vec(transpose(P))) * entr(dmtc(:, iOutput));
% 		for iInput = 1 : size(indexCombination, 2)
% 			combination = num2cell(indexCombination(:, iInput));
% 			z = jointArray(combination{:});
% 		end
	end





	for iTag = 1 : nTags
		for iState = 1 : nStates
			marginalSet = find(indexCombination(iTag, :) == iState);
			for iInput = marginalSet
				marginalInformation(iTag, iState, iInput) = prod(combinationDistribution(setdiff(tagSet, iTag), iInput), 1) * informationFunction(iInput);
			end
		end
	end
	rateBound = 1;
end
