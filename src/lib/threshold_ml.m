function [threshold, dmtc, backscatterRate] = threshold_ml(equivalentDistribution, receivedPower, symbolRatio, tolerance)
	% Function:
	%	- obtain the maximum likelihood detection thresholds based on typical power levels
    %
    % Input:
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination combination
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
    %
    % Output:
	%	- threshold [1 * (nOutputs + 1)] : the ML thresholding values
	%	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
	%	- tolerance: minimum non-zero input probability
    %
    % Comment:
    %	- ML detection thresholds do not depend on tag input distribution
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 17

	% * Declare default tolerance
	arguments
		equivalentDistribution;
		receivedPower;
		symbolRatio;
		tolerance = 1e-6;
	end

	% * Get data
% 	nOutputs = sum(equivalentDistribution >= tolerance);
	nOutputs = length(equivalentDistribution);

	% * Initialization
	[receivedPower, outputIndex] = sort(receivedPower);

	% * Each ML threshold depends on the average received power of adjacent detection regions
	threshold(nOutputs + 1) = inf;
	for iOutput = 2 : nOutputs
		threshold(iOutput) = symbolRatio * (receivedPower(iOutput - 1) * receivedPower(iOutput)) / (receivedPower(iOutput - 1) - receivedPower(iOutput)) * log(receivedPower(iOutput - 1) / receivedPower(iOutput));
	end

	% * Construct DMTC and compute mutual information
	dmtc = channel_discretization(threshold, receivedPower, symbolRatio);
	dmtc = dmtc(:, outputIndex);
	backscatterRate = rate_backscatter(equivalentDistribution, dmtc);
end
