function [thresholdDomain] = domain_threshold(symbolRatio, nBins, receivePower, confidence)
	% Function:
    %	- design (high-resolution) quantization set that act as domain of decision thresholds
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- nBins: number of receive energy quantization bins
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- confidence: confidence level for edge hypotheses
    %
    % Output:
	%	- thresholdDomain [1 x (nBins + 1)]: boundaries of quantized energy bins as domain of decision thresholds
    %
    % Comment:
	%	- discretize receive energy into numerous bins for threshold design
	%	- decision region of each letter consists of adjoint energy bins, and thresholds are selected from quantization set
	%	- high-resolution quantization improves detection performance at the cost of (linear) computation complexity
	%	- select quantizer uniformly over critical region associated with confidence level of edge hypotheses to improve performance
	%	- replace infinity by critical value of confidence approaching 1 (as infinity is not supported in PGD beamforming design)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 18

	% * Set default bin number and confidence level
	arguments
		symbolRatio;
		nBins;
		receivePower;
		confidence = 1 - 1e-3;
	end

	% * Compute lower, upper, and infinite critical thresholds
	tlb = icdf('Gamma', 1 - confidence, symbolRatio, min(receivePower));
	tub = icdf('Gamma', confidence, symbolRatio, max(receivePower));
	tib = icdf('Gamma', 1 - eps, symbolRatio, max(receivePower));

	% * Design threshold domain
	thresholdDomain = [0, linspace(tlb, tub, nBins - 1), tib];
end
