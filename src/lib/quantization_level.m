function [quantizationLevel] = quantization_level(symbolRatio, receivePower, nBins, confidence)
	% Function:
    %	- determine (high-resolution) quantization level set to quantize receive energy into numerous bins
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- nBins: number of receive energy quantization bins
    %	- confidence: number of standard deviations from mean (i.e., parameter k in Chebyshev's inequality)
    %
    % Output:
	%	- quantizationLevel [1 x (nBins + 1)]: boundaries of quantized energy bins
    %
    % Comment:
	%	- quantized receive energy bins should be grouped to formulate decision regions
	%	- hence, decision thresholds are selected from quantization threshold set
	%	- quantization threshold set is constructed within empirical interval by Chebyshev's inequality
	%	- high-resolution quantization threshold set improves detection performance at the cost of (linear) computation complexity
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 18

	% * Set default bin number and confidence score
	arguments
		symbolRatio;
		receivePower;
		nBins = 2 ^ 8;
		confidence = 3;
	end

	% * Obtain empirical lower and upper energy bounds by Chebyshev's inequality
	lowerBound = symbolRatio * min(receivePower) - confidence * sqrt(symbolRatio * min(receivePower));
	upperBound = symbolRatio * max(receivePower) + confidence * sqrt(symbolRatio * max(receivePower));

	% * Ensure non-negative thresholds
	if lowerBound <= 0
		quantizationLevel = [linspace(0, upperBound, nBins), Inf];
	else
		quantizationLevel = [0, linspace(lowerBound, upperBound, nBins - 1), Inf];
	end
end
