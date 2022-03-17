function [inputDistribution, equivalentDistribution, weightedSumRate] = recovery_marginalization(jointDistribution, dmtc, weight, symbolRatio, snr)
	% Function:
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by marginalization
    %
    % Input:
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight: the relative priority of the primary link
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
    %	- the simplest method to extract input distribution vectors from the joint input distribution array
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 09

	nTags = ndims_modified(jointDistribution);
	nStates = nthroot(size(dmtc, 1), nTags);
	inputDistribution = zeros(nTags, nStates);
	for iTag = 1 : nTags
		inputDistribution(iTag, :) = reshape(sum_nested(jointDistribution, setdiff(1 : nTags, iTag)), [1, nStates]);
	end
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc);
end
