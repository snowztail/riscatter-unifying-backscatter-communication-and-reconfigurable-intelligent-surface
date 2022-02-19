function [discreteChannel] = discretize_channel(threshold, receivedPower, symbolRatio)
	% Function:
	%	- obtain discrete channel based on channel probability distribution and bin boundaries
	%	- construct equivalent DMTC based on decision thresholds
    %
    % Input:
	%	- threshold [1 * nOutputs + 1]: bin boundaries or decision thresholds
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
    %
    % Output:
	%	- discreteChannel [(nStates ^ nTags) * nOutputs]: the DMC probability mass function after quantization (or the DMTC based on decision thresholds)
    %
    % Comment:
    %   - for a given tag input combination, the continous-output channel follows Erlang distribution
	%	- the discrete channel depends on bin boundaries or decision thresholds
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Get data
	nInputs = size(receivedPower, 1);
	nOutputs = size(threshold, 2) - 1;

	% * Construct DMC based on thresholding schemes
	discreteChannel = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			discreteChannel(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
end
