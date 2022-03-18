function [threshold, dmtc, backscatterRate] = threshold_dp(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio, tolerance)
	% Function:
	%	- obtain the DMTC capacity-achieving thresholding scheme by dynamic programming
    %
    % Input:
	%	- thresholdCandidate [1 * (nBins + 1)]: candidate threshold values
    %	- dmc [(nStates ^ nTags) * nBins]: the transition probability matrix of the backscatter discrete memoryless MAC obtained by quantization
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination combination
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
	%	- tolerance: minimum non-zero input probability
    %
    % Output:
	%	- threshold [1 * (nOutputs + 1)] : the optimal thresholding values
	%	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Declare default tolerance
	arguments
		thresholdCandidate;
		dmc;
		equivalentDistribution;
		receivedPower;
		symbolRatio;
		tolerance = eps;
	end

	% * Get data
% 	nOutputs = sum(equivalentDistribution >= tolerance);
	nOutputs = length(equivalentDistribution);
	nBins = size(dmc, 2);

	% * Initialization
	dp = zeros(nBins, nOutputs);
	sol = zeros(nBins, nOutputs);
	for iBin = 1 : nBins
		dp(iBin, 1) = quantization_cost(1 : iBin, equivalentDistribution, dmc);
	end

	% * Compute dp
	for iOutput = 2 : nOutputs
		for iBin = nBins - nOutputs + iOutput : - 1 : iOutput
			dpCandidate = inf(nBins - 2, 1);
			for iThreshold = iOutput - 1 : iBin - 1
				dpCandidate(iThreshold) = dp(iThreshold, iOutput - 1) + quantization_cost(iThreshold + 1 : iBin, equivalentDistribution, dmc);
			end
			[dp(iBin, iOutput), sol(iBin, iOutput)] = min(dpCandidate);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : - 1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdCandidate(index + 1);

	% * Construct DMTC and compute mutual information
	dmtc = channel_discretization(threshold, receivedPower, symbolRatio);
	backscatterRate = rate_backscatter(equivalentDistribution, dmtc);
end


function [quantizationCost] = quantization_cost(binIndex, equivalentDistribution, dmc)
	outputDistribution = equivalentDistribution * dmc;
	jointDistribution = transpose(equivalentDistribution) .* dmc;
	conditionalDistribution = sum(jointDistribution(:, binIndex), 2) / sum(outputDistribution(binIndex), 2);
	quantizationCost = - sum(outputDistribution(binIndex)) * sum(conditionalDistribution .* log(conditionalDistribution));
end
