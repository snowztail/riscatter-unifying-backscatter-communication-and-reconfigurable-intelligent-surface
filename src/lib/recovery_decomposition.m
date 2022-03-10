function [inputDistribution, equivalentDistribution, weightedSumRate] = recovery_decomposition(jointDistribution, dmtc, weight, symbolRatio, snr)
	% Function:
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by normalizing the best rank-1 CP approximation
    %
    % Input:
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
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
	weightedSumRate = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc);
end
