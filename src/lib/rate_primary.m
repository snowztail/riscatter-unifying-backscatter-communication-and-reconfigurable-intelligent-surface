function [primaryRate] = rate_primary(symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer)
	% (weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc)
	% Function:
    %	- compute the primary rate associated with each tag input combination
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %
    % Output:
	%	- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	% * Get data
	nInputs = size(equivalentDistribution, 2);

	% * Primary achievable rate
	snr = abs(equivalentChannel * beamformer) .^ 2 / noisePower;
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log(1 + snr(iInput));
	end
	primaryRate = equivalentDistribution * primaryInformationFunction;
end
