function [threshold] = thresholding(dmc, equivalentDistribution, thresholdCandidate)
	% Function:
	%	- obtain the DMTC capacity-achieving thresholding scheme by dynamic programming
    %
    % Input:
    %   - dmc [nInputs * nLevels]: the transition probability matrix of the backscatter discrete memoryless MAC obtained by quantization
	%	- equivalentDistribution [1 * nInputs]: optimal input combination probability distribution
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

	% * Compute dp
	for iOutput = 2 : nOutputs
		for iLevel = nLevels - nOutputs + iOutput : - 1 : iOutput
			dpCandidate = Inf(nLevels - 2, 1);
			for iThreshold = iOutput - 1 : iLevel - 1
				dpCandidate(iThreshold) = dp(iThreshold, iOutput - 1) + quantization_cost(iThreshold + 1, iLevel, equivalentDistribution, dmc);
			end
			[dp(iLevel, iOutput), sol(iLevel, iOutput)] = min(dpCandidate);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nLevels;
	for iOutput = nOutputs : - 1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdCandidate(index + 1);
end


function [cost] = quantization_cost(lowerLevel, upperLevel, equivalentDistribution, dmc)
	outputDistribution = equivalentDistribution * dmc;
	jointDistribution = transpose(equivalentDistribution) .* dmc;
	conditionalDistribution = sum(jointDistribution(:, lowerLevel : upperLevel), 2) / sum(outputDistribution(lowerLevel : upperLevel), 2);
	cost = - sum(outputDistribution(lowerLevel : upperLevel)) * sum(conditionalDistribution .* log2(conditionalDistribution));
end
