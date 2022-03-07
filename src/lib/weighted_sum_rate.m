function [weightedSumRate, primaryRate, backscatterRate] = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistribution, dmtc)
	% Function:
    %   - compute the weighted sum rate for a given input combination distribution and thresholding scheme
    %
    % Input:
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%	- primaryRate: the achievable rate for the primary link (bps/Hz)
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	primaryRate = equivalentDistribution * information_function_primary(symbolRatio, snr);
	backscatterRate = equivalentDistribution * information_function_backscatter(equivalentDistribution, dmtc);
	weightedSumRate = [primaryRate, backscatterRate] * weight;
end
