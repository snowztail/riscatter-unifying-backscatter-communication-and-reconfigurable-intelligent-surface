function [informationFunction] = information_function(weight, symbolRatio, snr, equivalentDistribution, dmtc)
	% Function:
	%	- compute the information function associated with each tag input combination
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- informationFunction [(nStates ^ nTags) * 1]: weighted sum information function associated with each input combination
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 19

	% * Get data
	[nInputs, nOutputs] = size(dmtc);

	% * Primary information function
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log(1 + snr(iInput));
	end

	% * Backscatter information function
	backscatterInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	backscatterInformationFunction = sum(backscatterInformationFunction, 2);

	% * Weighted sum information function
	informationFunction = weight * primaryInformationFunction + (1 - weight) * backscatterInformationFunction;
end
