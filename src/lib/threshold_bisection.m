function [threshold] = threshold_bisection(symbolRatio, receivePower, equivalentDistribution, thresholdDomain, binDmc, tolerance)
	% Function:
	%	- group received energy bins into convex adjoint decision regions by bisection
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- thresholdDomain [1 x (nBins + 1)]: boundaries of quantized energy bins as domain of decision thresholds	%	- binDmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
    %	- tolerance: maximum tolerable threshold/divergence error, and tail probability (i.e., 1 - confidence level) of replacing infinity by critical threshold
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
		thresholdDomain;
		binDmc;
		tolerance = eps;
	end

	% * Get data
	[nInputs, nBins] = size(binDmc);
	nOutputs = nInputs;

	% ! Fall back to initializer when underflow (backscatter rate approaching 0)
	persistent Initializer

	% * Optimal t_0 = 0 and t_{K+1} = ∞, traverse all possible t_1
	rate = 0;
	for iBin = 2 : nBins
		thresholdBisection = zeros(1, nOutputs + 1);
		thresholdBisection(2) = thresholdDomain(iBin);

		% * Update t_j for j = 2, ..., K - 1 consequently
		isValid = true;
		for iThreshold = 3 : nOutputs
			% * Reference threshold pair and divergence
			rtp = [thresholdBisection(iThreshold - 2), thresholdBisection(iThreshold - 1)];
			referenceDivergence = divergence_backward(symbolRatio, receivePower, equivalentDistribution, rtp);

			% * Lower and upper threshold pairs and divergences
			ltp = [thresholdBisection(iThreshold - 1), thresholdBisection(iThreshold - 1) + tolerance];
			utp = [thresholdBisection(iThreshold - 1), thresholdDomain(end)];
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
				if abs(middleDivergence - referenceDivergence) <= tolerance || (utp(end) - ltp(end)) <= tolerance
					thresholdBisection(iThreshold) = mtp(end);
					break;
				end
				if sign(middleDivergence - referenceDivergence) == sign(lowerDivergence - referenceDivergence)
					ltp = mtp;
					lowerDivergence = middleDivergence;
				else
					utp = mtp;
					% upperDivergence = middleDivergence;
				end
			end
		end
		% * Replace infinity by critical threshold of confidence approaching 1
		thresholdBisection(nOutputs + 1) = icdf('Gamma', 1 - tolerance, symbolRatio, max(receivePower));

		if isValid
			% * Valid threshold
			dmacBisection = dmc_integration(symbolRatio, receivePower, thresholdBisection);
			rateBisection = sum(entr(equivalentDistribution' * dmacBisection) - equivalentDistribution' * entr(dmacBisection));
			if rateBisection >= rate
				threshold = thresholdBisection;
				rate = rateBisection;
			end
		end
	end

	if exist('threshold', 'var')
		% * Store for fallback
		Initializer.threshold = threshold;
	else
		% * Invalid threshold (divergence approching 0), use previous result
		threshold = Initializer.threshold;
	end
end


function [divergence] = divergence_backward(symbolRatio, receivePower, equivalentDistribution, thresholdPair)
	% * Forward and backward bin probability
	fbp = dmc_integration(symbolRatio, receivePower, thresholdPair);
	bbp = equivalentDistribution .* fbp / (equivalentDistribution' * fbp);
	bbp(bbp < eps) = eps;

	% * Forward and backward threshold probability
	ftp = (thresholdPair(end) .^ (symbolRatio - 1) .* exp(-thresholdPair(end) ./ receivePower)) ./ (receivePower .^ symbolRatio .* gamma(symbolRatio));
	btp = equivalentDistribution .* ftp / (equivalentDistribution' * ftp);
	btp(btp < eps) = eps;

	% * Divergence of backward probabilities
	divergence = sum(rel_entr(btp, bbp));
end
