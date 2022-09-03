function [distribution, equivalentDistribution] = recovery_marginalization(jointDistribution)
	% Function:
	%	- extract individual tag input distribution from joint input distribution by marginalization
    %
    % Input:
	%	- jointDistribution [nStates x ... (nTags-dimensional) ... x nStates]: joint tag input distribution with full transmit cooperation
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 09

	% * Get data
	nTags = dimension_number(jointDistribution);
	nStates = size(jointDistribution, 1);

	% * Obtain marginal distribution
	distribution = zeros(nStates, nTags);
	for iTag = 1 : nTags
		distribution(:, iTag) = vec(sum(jointDistribution, [1 : iTag - 1, iTag + 1 : nTags]));
	end
	equivalentDistribution = prod(tuple_tag(distribution), 2);
end
