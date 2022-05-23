function [thresholdDomain] = domain_threshold(symbolRatio, receivePower, confidence, nBins)
	% Function:
    %	- determine (high-resolution) quantization levels to quantize receive energy into numerous bins
	%	- quantization levels act as domain of decision thresholds
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- confidence: coefficient to balance precision and underflow
	%	- nBins: number of receive energy quantization bins
    %
    % Output:
	%	- thresholdDomain [1 x (nBins + 1)]: boundaries of quantized energy bins as domain of decision thresholds
    %
    % Comment:
	%	- quantized receive energy bins should be grouped to formulate decision regions
	%	- hence, decision thresholds are selected from quantization levels
	%	- high-resolution quantization levels improves detection performance at the cost of (linear) computation complexity
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 18

	% * Set default bin number and confidence score
	arguments
		symbolRatio;
		receivePower;
		confidence = 10;
		nBins = 2 ^ 8;
	end

	% * Use ML threshold as reference
	threshold = threshold_ml(symbolRatio, receivePower);
	thresholdDomain = linspace(threshold(2) / confidence, threshold(end - 1) * confidence, nBins + 1);
end
