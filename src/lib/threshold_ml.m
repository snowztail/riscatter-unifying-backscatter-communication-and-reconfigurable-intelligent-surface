function [threshold] = threshold_ml(symbolRatio, receivePower)
	% Function:
	%	- obtain the decision thresholds for maximum likelihood detector
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of ML decision regions (including 0 and Inf)
    %
    % Comment:
    %	- ML detector does not depend on input distribution
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 17

	% * Get data
	nOutputs = size(receivePower, 1);

	% * Each ML threshold depends on average receive power of adjacent letters
	sortPower = sort(receivePower);
	threshold(nOutputs + 1) = Inf;
	for iOutput = 2 : nOutputs
		threshold(iOutput) = symbolRatio * (sortPower(iOutput - 1) * sortPower(iOutput)) / (sortPower(iOutput - 1) - sortPower(iOutput)) * log(sortPower(iOutput - 1) / sortPower(iOutput));
	end
end
