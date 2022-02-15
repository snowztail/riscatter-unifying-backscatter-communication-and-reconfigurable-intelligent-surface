function [threshold, sol] = thresholding_smawk(dmc, equivalentDistribution, thresholdCandidate)
	% Function:
	%	- obtain the DMTC capacity-achieving thresholding scheme by dynamic programming
    %
    % Input:
    %   - dmc [nInputs * nLevels]: the transition probability matrix of the backscatter discrete memoryless MAC obtained by quantization
	%	- equivalentDistribution: optimal input combination probability distribution
	%	- thresholdCandidate [1 * nLevels + 1]: all possible threshold values
    %
    % Output:
	%	- threshold: the optimal thresholding values
    %
    % Comment:
    %   - the dynamic programming algorithm can be accelerated by SWAWK once the quantization cost function satisfies the Quadrangle Inequality (QI)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Set default tolerance
	arguments
		dmc;
		equivalentDistribution;
		thresholdCandidate;
	end

	% * Ensure non-zero transitional probability as required by Blahut-Arimoto algorithm
	dmc(dmc < eps) = eps;

	% * Get data
	[nInputs, nLevels] = size(dmc);
	nOutputs = nInputs;

	% * Initialization
	dp = zeros(nLevels, nOutputs);
	sol = zeros(nLevels, nOutputs);
	for iLevel = 1 : nLevels
		dp(iLevel, 1) = quantization_cost(1, iLevel, equivalentDistribution, dmc);
	end

	% * Generate monotone matrix
	for iOutput = 2 : nOutputs
		D = zeros(nLevels - nOutputs + 1);
		for iRow = 1 : nLevels - nOutputs + 1
			for iColumn = 1 : nLevels - nOutputs + 1
				if iRow >= iColumn
					D(iRow, iColumn) = dp(iColumn - 2 + iOutput, iOutput - 1) + quantization_cost(iColumn - 1 + iOutput, iRow - 1 + iOutput, equivalentDistribution, dmc);
				else
					D(iRow, iColumn) = Inf;
				end
			end
		end

		% * Retrieve leftmost minima position by SMAWK
		[r, c] = deal(1 : nLevels - nOutputs + 1);
		[p] = smawk(D, r, c);

		% * Get sol and dp
		for iLevel = nLevels - nOutputs + iOutput : - 1 : iOutput
			sol(iLevel, iOutput) = p(iLevel - iOutput + 1) - 2 + iOutput;
			dp(iLevel, iOutput) = dp(sol(iLevel, iOutput), iOutput - 1) + quantization_cost(sol(iLevel, iOutput) + 1, iLevel, equivalentDistribution, dmc);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nLevels;
	for iOutput = nOutputs : - 1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdCandidate(index + 1);
end


function [quantizationCost] = quantization_cost(lowerLevel, upperLevel, equivalentDistribution, dmc)
	outputDistribution = equivalentDistribution * dmc;
	jointDistribution = transpose(equivalentDistribution) .* dmc;
	conditionalDistribution = sum(jointDistribution(:, lowerLevel : upperLevel), 2) / sum(outputDistribution(lowerLevel : upperLevel), 2);
	quantizationCost = - sum(outputDistribution(lowerLevel : upperLevel)) * sum(conditionalDistribution .* log2(conditionalDistribution));
end
