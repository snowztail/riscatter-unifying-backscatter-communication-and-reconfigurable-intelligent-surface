function [distribution, equivalentDistribution] = distribution_uniform(nTags, nStates)
	% Function:
	%	- obtain uniform input probability distribution
    %
    % Input:
	%	- nTags: number of tags
	%	- nStates: number of candidate reflecting states
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Comment:
    %	- used in uncoded backscatter transmission
    %
    % Author & Date: Yang (i@snowztail.com), 22 Sep 03

	distribution = normalize(ones(nStates, nTags), 'norm', 1);
	equivalentDistribution = prod(tuple_tag(distribution), 2);
end
