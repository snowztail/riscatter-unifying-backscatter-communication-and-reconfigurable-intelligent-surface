function [threshold, dmtc, backscatterRate] = t_smawk(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer)
	% Function:
	%	- group the received energy bins into convex decision regions by dynamic programming accelerated by SMAWK algorithm
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- dmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
	%	- thresholdSet [1 x (nBins + 1)]: boundaries of quantized energy bins
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
    % Comment:
    %	- SWAWK algorithm requires quantization cost to satisfy the Quadrangle Inequality (QI)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09


	% * Get data
	nOutputs = size(equivalentDistribution, 1);

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
					D(iRow, iColumn) = Inf;
				end
			end
		end

		% * Retrieve leftmost minima position by SMAWK
		[r, c] = deal(1 : nBins - nOutputs + 1);
		[p] = smawk(D, r, c);

		% * Get sol and dp
		for iBin = nBins - nOutputs + iOutput : -1 : iOutput
			sol(iBin, iOutput) = p(iBin - iOutput + 1) - 2 + iOutput;
			dp(iBin, iOutput) = dp(sol(iBin, iOutput), iOutput - 1) + quantization_cost(sol(iBin, iOutput) + 1 : iBin, equivalentDistribution, dmc);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : -1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdCandidate(index + 1);

	% TODO sort input to match threshold and DMTC?
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
