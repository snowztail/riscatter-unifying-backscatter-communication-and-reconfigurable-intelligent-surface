function [threshold, dmtc, backscatterRate] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio, tolerance)
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
	%   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
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
	backscatterRateSet = zeros(nLevels, 1);
	thresholdSet = cell(nLevels, 1);
	for iLevel = 2 : nLevels - 1
		threshold = zeros(1, nOutputs + 1);
		threshold(end) = inf;
		threshold(2) = thresholdCandidate(iLevel);

		% * Update t_j for j = 2, ..., K - 1 consequently
		isValid = true;
		for iThreshold = 3 : nOutputs
			% * Compute reference divergence
			binIndex = find(thresholdCandidate == threshold(iThreshold - 2)) : find(thresholdCandidate == threshold(iThreshold - 1)) - 1;
			thresholdDensity = threshold_density(receivedPower, symbolRatio, threshold(iThreshold - 1));
			divergence = backward_divergence(dmc, equivalentDistribution, binIndex, thresholdDensity);

			% * Initialize upper and lower bound of threshold
			lowerThreshold = thresholdCandidate(find(thresholdCandidate == threshold(iThreshold - 1)) + 1);
			upperThreshold = thresholdCandidate(end - 1);

			lowerBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == lowerThreshold) - 1;
			lowerDivergence = backward_divergence(dmc, equivalentDistribution, lowerBinIndex, thresholdDensity);
			lowerDeviation = lowerDivergence - divergence;

			upperBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == upperThreshold) - 1;
			upperDivergence = backward_divergence(dmc, equivalentDistribution, upperBinIndex, thresholdDensity);
			upperDeviation = upperDivergence - divergence;

			if (lowerDeviation < 0 && upperDeviation < 0) || (lowerDeviation > 0 && upperDeviation > 0) || isnan(divergence)
				isValid = false;
				threshold(iThreshold : end - 1) = nan;
				break;
			end

			% * Bisection to guarantee (approximately) equal divergence
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

				% * Shrink search interval and update deviations
				midThreshold = thresholdCandidate(round(mean([lowerBinIndex(end), upperBinIndex(end)])) + 1);
				midBinIndex = find(thresholdCandidate == threshold(iThreshold - 1)) : find(thresholdCandidate == midThreshold) - 1;
				midDivergence = backward_divergence(dmc, equivalentDistribution, midBinIndex, thresholdDensity);
				midDeviation = midDivergence - divergence;
				if midDeviation == 0
					threshold(iThreshold) = midThreshold;
					break;
				elseif (lowerDeviation < 0 && midDeviation < 0 && upperDeviation > 0) || (lowerDeviation > 0 && midDeviation > 0 && upperDeviation < 0)
					lowerThreshold = midThreshold;
					lowerBinIndex = midBinIndex;
					lowerDeviation = midDeviation;
				elseif (lowerDeviation < 0 && midDeviation > 0 && upperDeviation > 0) || (lowerDeviation > 0 && midDeviation < 0 && upperDeviation < 0)
					upperThreshold = midThreshold;
					upperBinIndex = midBinIndex;
					upperDeviation = midDeviation;
				end
			end
		end

		% * Construct DMTC and compute mutual information
		if isValid
			dmtc = discretize_channel(threshold, receivedPower, symbolRatio);
			backscatterRateSet(iLevel) = equivalentDistribution * information_function_backscatter(equivalentDistribution, dmtc);
			thresholdSet{iLevel} = threshold;
		else
			backscatterRateSet(iLevel) = nan;
			thresholdSet{iLevel} = threshold;
		end
	end
	[backscatterRate, optimalLevel] = max(backscatterRateSet);
	threshold = thresholdSet{optimalLevel};
	dmtc = discretize_channel(threshold, receivedPower, symbolRatio);
end


function [divergence] = backward_divergence(dmc, equivalentDistribution, binIndex, thresholdDensity)
	nInputs = size(equivalentDistribution, 2);
	backwardThreshold = zeros(nInputs, 1);
	backwardBin = zeros(nInputs, 1);
	divergence = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		backwardThreshold(iInput) = equivalentDistribution(iInput) * thresholdDensity(iInput) / (equivalentDistribution * thresholdDensity);
		backwardBin(iInput) = equivalentDistribution(iInput) * sum(dmc(iInput, binIndex), 2) / (equivalentDistribution * sum(dmc(:, binIndex), 2));
		divergence(iInput) = backwardThreshold(iInput) * log(backwardThreshold(iInput) / backwardBin(iInput));
	end
	divergence = sum(divergence, 1);
end

function [thresholdDensity] = threshold_density(receivedPower, symbolRatio, threshold)
	thresholdDensity = (threshold .^ (symbolRatio - 1) .* exp(-threshold ./ receivedPower)) ./ (receivedPower .^ symbolRatio .* gamma(symbolRatio));
end
