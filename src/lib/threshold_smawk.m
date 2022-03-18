function [threshold, dmtc, backscatterRate] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio, tolerance)
	% Function:
	%	- obtain the DMTC capacity-achieving thresholding scheme by dynamic programming accelerated by SMAWK algorithm
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
    % Comment:
    %	- SWAWK algorithm requires the quantization cost function to satisfy the Quadrangle Inequality (QI)
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
	dmtc = channel_discretization(threshold, receivedPower, symbolRatio);
	backscatterRate = rate_backscatter(equivalentDistribution, dmtc);
end


function [quantizationCost] = quantization_cost(binIndex, equivalentDistribution, dmc)
	outputDistribution = equivalentDistribution * dmc;
	jointDistribution = transpose(equivalentDistribution) .* dmc;
	conditionalDistribution = sum(jointDistribution(:, binIndex), 2) / sum(outputDistribution(binIndex), 2);
	quantizationCost = - sum(outputDistribution(binIndex)) * sum(conditionalDistribution .* log(conditionalDistribution));
end
