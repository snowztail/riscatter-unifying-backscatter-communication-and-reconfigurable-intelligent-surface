function [threshold] = threshold_ml(symbolRatio, receivePower, tolerance)
	% Function:
	%	- obtain the decision thresholds for maximum likelihood detector
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- tolerance: tail probability (i.e., 1 - confidence level) of replacing infinity by critical threshold
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of ML decision regions (including 0 and Inf)
    %
    % Comment:
	%	- for completness, zero and infinity are defined as first and last thresholds
    %	- other thresholds depend on average receiver power of (power-sorted) adjacent letters, and
	%	- threshold design is not influenced by input distribution
	%	- replace infinity by critical threshold of confidence approaching 1 (as infinity is not supported in PGD beamforming design)
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 17

	% * Set default tolerance
	arguments
		symbolRatio;
		receivePower;
		tolerance = eps;
	end

	% * Get data
	nOutputs = size(receivePower, 1);

	% * Design ML thresholds
	sortPower = sort(receivePower);
	% threshold(1) = 0;
	threshold(nOutputs + 1) = icdf('Gamma', 1 - tolerance, symbolRatio, max(receivePower));
	for iOutput = 2 : nOutputs
		threshold(iOutput) = symbolRatio * (sortPower(iOutput - 1) * sortPower(iOutput)) / (sortPower(iOutput - 1) - sortPower(iOutput)) * log(sortPower(iOutput - 1) / sortPower(iOutput));
	end
end
