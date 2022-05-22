function [stateTuple] = state_tuple(nTags, nStates)
	% Function:
	%	- enumerate tag state tuples with entries as state indexes
    %
    % Input:
	%	- nTags: number of tags
    %	- nStates: number of available states at tags (i.e., modulation order)
    %
    % Output:
	%	- stateTuple [nTags x nInputs]: each column is a nTags-tuple whose entries come from state indexes
    %
    % Author & Date: Yang (i@snowztail.com), 21 Dec 25

	% * Enumerate tuples by recursion
	if nTags > 1
		stateTuple = sortrows(combvec(1 : nStates, state_tuple(nTags - 1, nStates)));
	elseif nTags == 1
		stateTuple = 1 : nStates;
	% elseif nTags == 0
	% 	stateTuple = double.empty(0, 1);
	end
end
