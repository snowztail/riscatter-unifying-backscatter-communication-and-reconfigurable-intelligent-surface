function [inputDistribution, equivalentDistribution, weightedSumRate] = recovery_marginalization(weight, symbolRatio, equivalentChannel, noisePower, jointDistribution, beamformer, dmtc)
	% Function:
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by marginalization
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
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
    %	- obtain the marginal distribution from the joint distribution array
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 09

	% * Get data
	nTags = ndims_modified(jointDistribution);
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Obtain marginal distribution
	inputDistribution = zeros(nTags, nStates);
	for iTag = 1 : nTags
		inputDistribution(iTag, :) = reshape(sum_nested(jointDistribution, setdiff(1 : nTags, iTag)), [1, nStates]);
	end
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
end
