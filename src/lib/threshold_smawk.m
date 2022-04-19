function [threshold, dmtc, backscatterRate] = threshold_smawk(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer)
	% Function:
	%	- group the received energy bins into convex decision regions by dynamic programming accelerated by SMAWK algorithm
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
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
    % Comment:
    %	- SWAWK algorithm requires the quantization cost function to satisfy the Quadrangle Inequality (QI)
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

	% * Generate monotone matrix
	for iOutput = 2 : nOutputs
		D = zeros(nBins - nOutputs + 1);
		for iRow = 1 : nBins - nOutputs + 1
			for iColumn = 1 : nBins - nOutputs + 1
				if iRow >= iColumn
					D(iRow, iColumn) = dp(iColumn - 2 + iOutput, iOutput - 1) + quantization_cost(iColumn - 1 + iOutput : iRow - 1 + iOutput, equivalentDistribution, dmc);
				else
					D(iRow, iColumn) = inf;
				end
			end
		end

		% * Retrieve leftmost minima position by SMAWK
		[r, c] = deal(1 : nBins - nOutputs + 1);
		[p] = smawk(D, r, c);

		% * Get sol and dp
		for iBin = nBins - nOutputs + iOutput : - 1 : iOutput
			sol(iBin, iOutput) = p(iBin - iOutput + 1) - 2 + iOutput;
			dp(iBin, iOutput) = dp(sol(iBin, iOutput), iOutput - 1) + quantization_cost(sol(iBin, iOutput) + 1 : iBin, equivalentDistribution, dmc);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : - 1 : 1
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
