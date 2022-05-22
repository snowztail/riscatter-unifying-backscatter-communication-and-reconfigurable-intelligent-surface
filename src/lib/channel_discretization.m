function [discreteChannel, receivedPower] = channel_discretization(symbolRatio, equivalentChannel, noisePower, beamformer, threshold)
	% Function:
	%	- obtain energy detection channel based on channel probability distribution and bin boundaries
	%	- construct equivalent DMTC based on decision thresholds
    %
    % Input:
	%	- symbolRatio: the ratio of the secondary symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent primary channel under each tag input combination
	%	- noisePower: average noise power at the user
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
	%	- threshold [1 * nOutputs + 1]: bin boundaries or decision thresholds
    %
    % Output:
	%	- discreteChannel [(nStates ^ nTags) * nOutputs]: the DMC probability mass function after quantization (or the DMTC based on decision thresholds)
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination combination
    %
    % Comment:
    %	- for a given tag input combination, the continous output channel follows Erlang distribution
	%	- the discrete channel depends on bin boundaries or decision thresholds
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Get data
	nInputs = size(equivalentChannel, 1);
	nOutputs = size(threshold, 2) - 1;

	% * Compute the expected received power under all backscatter input combinations
	if isvector(beamformer)
		% * Beamformer vector
		receivedPower = abs(equivalentChannel * beamformer) .^ 2 + noisePower;
	else
		% * Beamformer matrix
		receivedPower = zeros(nInputs, 1);
		for iInput = 1 : nInputs
			receivedPower(iInput) = real(trace(equivalentChannel(iInput, :)' * equivalentChannel(iInput, :) * beamformer) + noisePower);
		end
	end

	% * Construct DMC based on thresholding schemes
	discreteChannel = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			discreteChannel(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
end
