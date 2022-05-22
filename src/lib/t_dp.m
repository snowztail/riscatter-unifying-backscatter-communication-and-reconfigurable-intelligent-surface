function [threshold, dmtc, backscatterRate] = threshold_dp(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer)
	% Function:
	%	- group the received energy bins into convex decision regions by dynamic programming
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent primary channel under each tag input combination
	%	- noisePower: average noise power at the user
	%	- nBins: number of discretization bins over received signal
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %
    % Output:
	%	- threshold [1 * (nOutputs + 1)]: boundaries of decision regions
	%	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Get data
	nOutputs = size(equivalentDistribution, 2);

	% * Obtain threshold candidates that delimit output into discrete bins
	thresholdCandidate = threshold_candidate(symbolRatio, equivalentChannel, noisePower, nBins, beamformer);

	% * Evaluate DMC over all bins
	dmc = channel_discretization(symbolRatio, equivalentChannel, noisePower, beamformer, thresholdCandidate);

	% * Initialize cost functions
	dp = zeros(nBins, nOutputs);
	sol = zeros(nBins, nOutputs);
	for iBin = 1 : nBins
		dp(iBin, 1) = quantization_cost(1 : iBin, equivalentDistribution, dmc);
	end

	% * Compute dp
	for iOutput = 2 : nOutputs
		for iBin = nBins - nOutputs + iOutput : -1 : iOutput
			dpCandidate = inf(nBins - 2, 1);
			for iThreshold = iOutput - 1 : iBin - 1
				dpCandidate(iThreshold) = dp(iThreshold, iOutput - 1) + quantization_cost(iThreshold + 1 : iBin, equivalentDistribution, dmc);
			end
			[dp(iBin, iOutput), sol(iBin, iOutput)] = min(dpCandidate);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : -1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdCandidate(index + 1);

	% * Construct DMTC and compute mutual information
	dmtc = channel_discretization(symbolRatio, equivalentChannel, noisePower, beamformer, threshold);
	backscatterRate = rate_backscatter(equivalentDistribution, dmtc);
end


function [quantizationCost] = quantization_cost(binIndex, equivalentDistribution, dmc)
	outputDistribution = equivalentDistribution * dmc;
	jointDistribution = transpose(equivalentDistribution) .* dmc;
	conditionalDistribution = sum(jointDistribution(:, binIndex), 2) / sum(outputDistribution(binIndex), 2);
	quantizationCost = - sum(outputDistribution(binIndex)) * sum(conditionalDistribution .* log(conditionalDistribution));
end
