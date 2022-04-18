function [inputDistribution, equivalentDistribution, weightedSumRate] = recovery_decomposition(weight, symbolRatio, equivalentChannel, noisePower, jointDistribution, precoder, dmtc)
	% Function:
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by normalizing the best rank-1 CP approximation
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
	%	- precoder [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
    %	- approximate the joint input distribution tensor by the best rank-1 CP tensor
	%
	% Reference:
	%	- General software, latest release: Brett W. Bader, Tamara G. Kolda and others, Tensor Toolbox for MATLAB, Version 3.2.1, www.tensortoolbox.org, April 5, 2021.
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 10

	cpTensor = cp_als(tensor(jointDistribution), 1);
	inputDistribution = cell2mat(cellfun(@transpose, cpTensor.U, 'UniformOutput', false));
	inputDistribution = inputDistribution ./ sum(inputDistribution, 2);
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, precoder, dmtc);
end
