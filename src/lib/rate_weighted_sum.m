function [weightedSumRate, primaryRate, backscatterRate] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc)
	% Function:
    %	- compute the weighted sum rate for a given input combination distribution and thresholding scheme
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%	- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22


	primaryRate = rate_primary(symbolRatio, snr, equivalentDistribution);
	backscatterRate = rate_backscatter(equivalentDistribution, dmtc);
	weightedSumRate = weight * primaryRate + (1 - weight) * backscatterRate;
end