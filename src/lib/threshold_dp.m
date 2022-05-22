function [threshold] = threshold_dp(equivalentDistribution, quantizationSet, binDmc)
	% Function:
	%	- group received energy bins into convex adjoint decision regions by dynamic programming
    %
    % Input:
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- quantizationSet [1 x (nBins + 1)]: boundaries of quantized energy bins
	%	- binDmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
	% Reference:
	%	- X. He, K. Cai, W. Song, and Z. Mei, “Dynamic programming for sequential deterministic quantization of discrete memoryless channels,” IEEE Transactions on Communications, vol. 69, no. 6, pp. 3638–3651, 2021.
	%
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Get data
	[nInputs, nBins] = size(binDmc);
	nOutputs = nInputs;

	% * Obtain relevant distributions
	jointDistribution = equivalentDistribution .* binDmc;
	outputDistribution = equivalentDistribution' * binDmc;

	% * Initialize cost functions
	dp = zeros(nBins, nOutputs);
	sol = zeros(nBins, nOutputs);
	for iBin = 1 : nBins
		dp(iBin, 1) = cost_quantization(1 : iBin, jointDistribution, outputDistribution);
	end

	% * Compute dp
	for iOutput = 2 : nOutputs
		for iBin = nBins - nOutputs + iOutput : -1 : iOutput
			dpCandidate = inf(nBins - 2, 1);
			for iThreshold = iOutput - 1 : iBin - 1
				dpCandidate(iThreshold) = dp(iThreshold, iOutput - 1) + cost_quantization(iThreshold + 1 : iBin, jointDistribution, outputDistribution);
			end
			[dp(iBin, iOutput), sol(iBin, iOutput)] = min(dpCandidate);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : -1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = quantizationSet(index + 1);
end
