function [rateBound] = rate_bound(dmtc, subSet, jointArray, complementArray)
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

% 	jointArray = rand(size(jointArray));
% 	complementArray = rand(size(complementArray));
	% * Get data
	nTags = ndims(jointArray);
	nStates = size(jointArray, 1);
	nOutputs = size(dmtc, 2);

	% * Initialization
	setIndexCombination = combvec_nested(1 : nStates, nTags);
	complementSet = setdiff(1 : nTags, subSet);
	subSetIndexCombination = combvec_nested(1 : nStates, length(subSet));
	complementSetIndexCombination = combvec_nested(1 : nStates, nTags - length(subSet));
	nSubSets = size(subSetIndexCombination, 2);
	nComplementSets = size(complementSetIndexCombination, 2);

	% * Formulate upper bounds
	jointDistribution = cvx(zeros(nOutputs, nComplementSets, nSubSets));
	conditionalEntropy = cvx(zeros(nOutputs, 1));
	relativeEntropy = cvx(zeros(nOutputs, nComplementSets));
	for iOutput = 1 : nOutputs
		% * Simplify operation by rearranging array dimensions in natural order
		conditionalEntropy(iOutput) = transpose(vec(permute(jointArray, ndims(jointArray) : -1 : 1))) * entr(dmtc(:, iOutput));
		for iComplementSet = 1 : nComplementSets
			complementSetIndex = transpose(complementSetIndexCombination(:, iComplementSet));
			complementSetIndexCell = num2cell(complementSetIndex);
			for iSubSet = 1 : nSubSets
				subSetIndex = transpose(subSetIndexCombination(:, iSubSet));
				setIndex = [subSetIndex complementSetIndex];
				setIndex([subSet complementSet]) = setIndex;
				setIndexCell = num2cell(setIndex);
				jointDistribution(iOutput, iComplementSet, iSubSet) = jointArray(setIndexCell{:}) * dmtc(iOutput, all(setIndexCombination == transpose(setIndex)));
			end
			relativeEntropy(iOutput, iComplementSet) = rel_entr(sum(jointDistribution(iOutput, iComplementSet, :)), complementArray(complementSetIndexCell{:}));
		end
	end
	rateBound = - sum(conditionalEntropy) - sum(sum(relativeEntropy));


% 	for iOutput = 1 : nOutputs
% 		pp12 = transpose(vec(transpose(P))) * dmtc(:, iOutput);
% % 			f12t2(iOutput) = rel_entr(pp12, 1);
% 		f12t2(iOutput) = - entr(pp12);
% 		ft1(iOutput) = transpose(vec(transpose(P))) * entr(dmtc(:, iOutput));
% 
% 		for iState = 1 : nStates
% 			iIndex = indexCombination(1, :) == iState;
% 			pp2 = P(iState, :) * dmtc(iIndex, iOutput);
% 			f2t2(iState, iOutput) = rel_entr(pp2, p1(iState));
% 		end
% 
% 		for jState = 1 : nStates
% 			jIndex = indexCombination(2, :) == jState;
% 			pp1 = transpose(P(:, jState)) * dmtc(jIndex, iOutput);
% 			f1t2(jState, iOutput) = rel_entr(pp1, p2(jState));
% 		end
% 	end
% 	f1 = - sum(ft1) - sum(sum(f1t2));
% 	f2 = - sum(ft1) - sum(sum(f2t2));
% 	f12 = - sum(ft1) - sum(f12t2);
% 
% 
% 	for iTag = 1 : nTags
% 		for iState = 1 : nStates
% 			marginalSet = find(indexCombination(iTag, :) == iState);
% 			for iInput = marginalSet
% 				marginalInformation(iTag, iState, iInput) = prod(combinationDistribution(setdiff(tagSet, iTag), iInput), 1) * informationFunction(iInput);
% 			end
% 		end
% 	end
% 	rateBound = 1;
end
