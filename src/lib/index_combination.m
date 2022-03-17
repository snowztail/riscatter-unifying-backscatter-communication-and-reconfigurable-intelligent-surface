function [indexCombination] = index_combination(nTags, nStates)
	% Function:
    %	- create all combinations of copies of a indexVector
    %
    % Input:
	%	- nTags: number of tags
    %	- nStates: number of states in tag constellation diagram
    %
    % Output:
	%	- indexCombination [nTags * (nStates ^ nTags)]: each column represents a possible tag index combination
    %
    % Author & Date: Yang (i@snowztail.com), 21 Dec 25

	indexVector = 1 : nStates;
	if nTags > 1
		indexCombination = sortrows(combvec(indexVector, index_combination(nTags - 1, nStates)));
	elseif nTags == 1
		indexCombination = indexVector;
	elseif nTags == 0
		indexCombination = double.empty(0, 1);
	end
end
