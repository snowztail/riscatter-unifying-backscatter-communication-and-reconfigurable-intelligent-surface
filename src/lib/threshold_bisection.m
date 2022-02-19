function [threshold, backscatterMutualInformation] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio, tolerance)
	% Function:
	%	- obtain the DMTC capacity-achieving thresholding scheme by bisection
    %
    % Input:
	%	- thresholdCandidate [1 * (nBins + 1)]: candidate threshold values
    %   - dmc [(nStates ^ nTags) * nBins]: the transition probability matrix of the backscatter discrete memoryless MAC obtained by quantization
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination combination
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
	%	- tolerance: maximum bisection deviation from zero
    %
    % Output:
	%	- threshold [1 * nOutputs] : the optimal thresholding values
	%	- backscatterMutualInformation: the optimal sum backscatter mutual information
    %
    % Comment:
    %   - For a given t_0 and t_1, the optimal t_j for j = 2, ..., K - 1 can be uniquely determined by bisection
	%	- Traverse all possible t_1 and choose the corresponding thresholding set that maximizes mutual information
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 18

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
	nOutputs = sum(equivalentDistribution >= tolerance);
	nLevels = size(thresholdCandidate, 2);

	% * Optimal t_0 = 0 and t_K+1 = ∞, traverse all possible t_1
	backscatterMutualInformation = zero(nLevels, 1);
	thresholdSet = cell(nLevels, 1);
	for iLevel = 1 : nLevels
		threshold = zeros(1, nOutputs + 1);
		threshold(end) = Inf;
		threshold(2) = thresholdCandidate(iLevel);

		% * Update t_j for j = 2, ..., K - 1 consequently
		for iThreshold = 3 : nOutputs
			% * Compute reference divergence
			binIndex = find(thresholdCandidate == threshold(iThreshold - 2)) : find(thresholdCandidate == threshold(iThreshold - 1));
			binDensity = dmc(:, thresholdCandidate == threshold(iThreshold - 1));
			divergence = backward_divergence(dmc, equivalentDistribution, binIndex, binDensity);

			% * Bisection to guarantee (approximately) equal divergence
			lowerThreshold = threshold(iThreshold - 1);
			lowerBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == lowerThreshold);
			lowerDivergence = backward_divergence(dmc, equivalentDistribution, lowerBinIndex, binDensity);
			lowerDeviation = lowerDivergence - divergence;

			upperThreshold = thresholdCandidate(end);
			upperBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == upperThreshold);
			upperDivergence = backward_divergence(dmc, equivalentDistribution, upperBinIndex, binDensity);
			upperDeviation = upperDivergence - divergence;

			if (lowerDeviation < 0 && upperDeviation < 0) || (lowerDeviation > 0 && upperDeviation > 0)
				error('Invalid initialization.');
			end

			while true
				% * Tolerable deviation, exit loop and update threshold
				if abs(lowerDeviation) <= tolerance || abs(upperDeviation) <= tolerance || find(thresholdCandidate == upperThreshold) - find(thresholdCandidate == lowerThreshold) <= 1
					if abs(lowerDeviation) <= abs(upperDeviation)
						threshold(iThreshold) = lowerThreshold;
					else
						threshold(iThreshold) = upperThreshold;
					end
					break;
				end

				% * Shrink search interval and continue looping
				[~, midIndex] = min(abs(thresholdCandidate - mean([lowerThreshold, upperThreshold])));
				midThreshold = threshold(midIndex);
				midBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == midThreshold);
				midDivergence = backward_divergence(dmc, equivalentDistribution, midBinIndex, binDensity);
				midDeviation = midDivergence - divergence;
				if midDeviation == 0
					threshold(iThreshold) = midThreshold;
					break;
				elseif (lowerDeviation < 0 && midDeviation < 0 && upperDeviation > 0) || (lowerDeviation > 0 && midDeviation > 0 && upperDeviation < 0)
					lowerThreshold = midThreshold;
				elseif (lowerDeviation < 0 && midDeviation > 0 && upperDeviation > 0) || (lowerDeviation > 0 && midDeviation < 0 && upperDeviation < 0)
					upperThreshold = midThreshold;
				end

				% * Update deviations
				lowerBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == lowerThreshold);
				lowerDivergence = backward_divergence(dmc, equivalentDistribution, lowerBinIndex, binDensity);
				lowerDeviation = lowerDivergence - divergence;

				upperBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == upperThreshold);
				upperDivergence = backward_divergence(dmc, equivalentDistribution, upperBinIndex, binDensity);
				upperDeviation = upperDivergence - divergence;
			end
		end

		% * Construct DMTC and compute mutual information
		dmtc = discretize_channel(threshold, receivedPower, symbolRatio);
		backscatterMutualInformation(iLevel) = equivalentDistribution * information_function_backscatter(equivalentDistribution, dmtc);
		thresholdSet{iLevel} = threshold;
	end
	[backscatterMutualInformation, optimalLevel] = max(backscatterMutualInformation);
	threshold = thresholdSet{optimalLevel};
end


function [divergence] = backward_divergence(dmc, equivalentDistribution, binIndex, binDensity)
	nInputs = size(equivalentDistribution, 2);
	backwardBin = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		backwardBin(iInput) = equivalentDistribution(iInput) * binDensity(iInput) / (equivalentDistribution * binDensity);
	end
	backwardRegion = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		backwardRegion(iInput) = equivalentDistribution(iInput) * sum(dmc(iInput, binIndex), 2) / (equivalentDistribution * sum(dmc(:, binIndex), 2));
	end
	divergence = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		divergence(iInput) = backwardBin(iInput) * log(backwardBin(iInput) / backwardRegion(iInput));
	end
	divergence = sum(divergence, 1);
end
