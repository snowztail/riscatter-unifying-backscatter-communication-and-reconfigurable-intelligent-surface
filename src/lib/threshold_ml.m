function [threshold, dmtc] = threshold_ml(symbolRatio, equivalentChannel, noisePower, precoder)
	% Function:
	%	- obtain the decision thresholds of maximum likelihood detector
    %
    % Input:
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- precoder [nTxs * 1]: transmit beamforming vector at the AP
    %
    % Output:
	%	- threshold [1 * (nOutputs + 1)] : the ML thresholding values
	%	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Comment:
    %	- ML detection does not require input distribution
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 17

	% * Get data
	nOutputs = size(equivalentChannel, 1);

	% * Compute the expected received power under all backscatter input combinations
	[sortedPower, outputIndex] = sort(abs(equivalentChannel * precoder) .^ 2 + noisePower);

	% * Each ML threshold depends on the average received power of adjacent detection regions
	threshold(nOutputs + 1) = inf;
	for iOutput = 2 : nOutputs
		threshold(iOutput) = symbolRatio * (sortedPower(iOutput - 1) * sortedPower(iOutput)) / (sortedPower(iOutput - 1) - sortedPower(iOutput)) * log(sortedPower(iOutput - 1) / sortedPower(iOutput));
	end

	% * Construct DMTC
	dmtc = channel_discretization(symbolRatio, equivalentChannel, noisePower, precoder, threshold);
	dmtc = dmtc(:, outputIndex);
end
