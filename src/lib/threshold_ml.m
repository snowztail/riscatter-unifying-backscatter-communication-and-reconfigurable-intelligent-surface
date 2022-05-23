function [threshold] = threshold_ml(symbolRatio, receivePower, confidence)
	% Function:
	%	- obtain the decision thresholds for maximum likelihood detector
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- confidence: coefficient to balance precision and underflow
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of ML decision regions (including 0 and Inf)
    %
    % Comment:
    %	- ML detector does not depend on input distribution
	%	- replace infinity by finite threshold to avoid precision issue (underflow)
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 17

	% * Set default confidence score
	arguments
		symbolRatio;
		receivePower;
		confidence = 10;
	end

	% * Get data
	nOutputs = size(receivePower, 1);

	% * Each ML threshold depends on average receive power of adjacent letters
	sortPower = sort(receivePower);

	% ! Replace infinity by finite threshold to avoid precision issue (underflow)
	threshold = zeros(1, nOutputs + 1);
	for iOutput = 2 : nOutputs
		threshold(iOutput) = symbolRatio * (sortPower(iOutput - 1) * sortPower(iOutput)) / (sortPower(iOutput - 1) - sortPower(iOutput)) * log(sortPower(iOutput - 1) / sortPower(iOutput));
	end
	% threshold(nOutputs + 1) = Inf;
	threshold(end) = confidence * threshold(end - 1);
end
