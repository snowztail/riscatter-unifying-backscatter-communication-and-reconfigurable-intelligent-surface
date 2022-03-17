function [primaryRate] = rate_primary(symbolRatio, snr, equivalentDistribution)
	% Function:
    %	- compute the primary rate associated with each tag input combination
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %
    % Output:
	%	- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	nInputs = size(snr, 1);
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log(1 + snr(iInput));
	end
	primaryRate = equivalentDistribution * primaryInformationFunction;
end
