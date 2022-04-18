function [weightedSumRate, rate] = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, precoder, dmtc)
	% Function:
    %	- compute the weighted sum rate for a given input combination distribution and thresholding scheme
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- precoder [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%	- rate [2 * 1]:
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	rate = [rate_primary(symbolRatio, equivalentChannel, noisePower, equivalentDistribution, precoder); rate_backscatter(equivalentDistribution, dmtc)];
	weightedSumRate = [weight, 1 - weight] * rate;
end
