function [distribution, equivalentDistribution] = recovery_decomposition(jointDistribution)
	% Function:
	%	- extract individual tag input distribution from joint input distribution by normalizing the best rank-1 CP approximation
    %
    % Input:
	%	- jointDistribution [nStates x ... (nTags-dimensional) ... x nStates]: joint tag input distribution with full transmit cooperation
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Comment:
    %	- approximate the joint input distribution tensor by the best rank-1 CP tensor
	%
	% Reference:
	%	- General software, latest release: Brett W. Bader, Tamara G. Kolda and others, Tensor Toolbox for MATLAB, Version 3.2.1, www.tensortoolbox.org, April 5, 2021.
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 10

	cpTensor = cp_als(tensor(jointDistribution), 1);
	distribution = normalize(cell2mat(transpose(cpTensor.U)), 'norm', 1);
	equivalentDistribution = prod(tuple_tag(distribution), 2);
end
