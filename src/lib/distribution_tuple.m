function [distributionTuple] = distribution_tuple(distribution)
	% Function:
    %	- obtain all combinations of tag input probability
    %
    % Input:
	%	- distribution [nTags * nStates]: input probability distribution
    %
    % Output:
	%	- combinationDistribution [nTags * (nStates ^ nTags)]: each column represents a possible tag input probability combination
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	[nTags, nStates] = size(distribution);
	if nTags > 1
% 		distributionTuple = flipud(combvec(distribution(1, :), distribution_tuple(distribution(2 : end, :))));
		distributionTuple = combvec(distribution_tuple(distribution(1 : end - 1, :)), distribution(end, :));
	elseif nTags == 1
		distributionTuple = distribution;
	end
end
