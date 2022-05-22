function [threshold] = threshold_bisection(symbolRatio, receivePower, equivalentDistribution, quantizationSet, binDmc, tolerance)
	% Function:
	%	- group received energy bins into convex adjoint decision regions by bisection
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- quantizationSet [1 x (nBins + 1)]: boundaries of quantized energy bins
	%	- binDmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
    % Comment:
    %	- for a given t_0 and t_1, the optimal t_j for j = 2, ..., K - 1 can be uniquely determined by bisection
	%	- traverse all possible t_1 and choose the threshold set that maximizes mutual information
    %
	% Reference:
	%	- T. Nguyen and T. Nguyen, “On thresholding quantizer design for mutual information maximization: Optimal Structures and algorithms,” 2020 IEEE 91st Vehicular Technology Conference (VTC2020-Spring), 2020.
	%
    % Author & Date: Yang (i@snowztail.com), 22 Feb 18

	% * Set default tolerance
	arguments
		symbolRatio;
		receivePower;
		equivalentDistribution;
		quantizationSet;
		binDmc;
		tolerance = 1e-6;
	end

	% * Get data
	[nInputs, nBins] = size(binDmc);
	nOutputs = nInputs;

	% * Optimal t_0 = 0 and t_{K+1} = ∞, traverse all possible t_1
	rate = 0;
	for iBin = 2 : nBins
		thresholdBisection = zeros(1, nOutputs + 1);
		[thresholdBisection(2), thresholdBisection(nOutputs + 1)] = deal(quantizationSet(iBin), Inf);

		% * Update t_j for j = 2, ..., K - 1 consequently
		isValid = true;
		for iThreshold = 3 : nOutputs
			% * Reference threshold pair and divergence
			rtp = [thresholdBisection(iThreshold - 2), thresholdBisection(iThreshold - 1)];
			referenceDivergence = divergence_backward(symbolRatio, receivePower, equivalentDistribution, rtp);

			% * Lower and upper threshold pairs and divergences
			ltp = [thresholdBisection(iThreshold - 1), thresholdBisection(iThreshold - 1) + tolerance];
			utp = [thresholdBisection(iThreshold - 1), quantizationSet(end - 1)];
			lowerDivergence = divergence_backward(symbolRatio, receivePower, equivalentDistribution, ltp);
			upperDivergence = divergence_backward(symbolRatio, receivePower, equivalentDistribution, utp);

			% * Search next threshold by bisection
			if sign(lowerDivergence - referenceDivergence) == sign(upperDivergence - referenceDivergence)
				isValid = false;
				break;
			end
			% assert(sign(lowerDivergence - referenceDivergence) ~= sign(upperDivergence - referenceDivergence));
			while true
				mtp = 0.5 * (ltp + utp);
				middleDivergence = divergence_backward(symbolRatio, receivePower, equivalentDistribution, mtp);
				if abs(middleDivergence - referenceDivergence) <= tolerance || (ltp(end) - utp(end)) <= tolerance
					thresholdBisection(iThreshold) = mtp(end);
					break;
				end
				if sign(middleDivergence - referenceDivergence) == sign(lowerDivergence - referenceDivergence)
					ltp = mtp;
					lowerDivergence = middleDivergence;
				else
					utp = mtp;
					upperDivergence = middleDivergence;
				end
			end
		end

		% * Check optimality
		if isValid
			dmacBisection = dmc_integration(symbolRatio, receivePower, thresholdBisection);
			rateBisection = sum(entr(equivalentDistribution' * dmacBisection) - equivalentDistribution' * entr(dmacBisection));
			if rateBisection >= rate
				threshold = thresholdBisection;
				rate = rateBisection;
			end
		end
	end
end


function [divergence] = divergence_backward(symbolRatio, receivePower, equivalentDistribution, thresholdPair)
	% * Forward and backward bin probability
	fbp = dmc_integration(symbolRatio, receivePower, thresholdPair);
	bbp = equivalentDistribution .* fbp / (equivalentDistribution' * fbp);

	% * Forward and backward threshold probability
	ftp = (thresholdPair(end) .^ (symbolRatio - 1) .* exp(-thresholdPair(end) ./ receivePower)) ./ (receivePower .^ symbolRatio .* gamma(symbolRatio));
	btp = equivalentDistribution .* ftp / (equivalentDistribution' * ftp);

	% * Divergence of backward probabilities
	divergence = sum(rel_entr(btp, bbp));
end
