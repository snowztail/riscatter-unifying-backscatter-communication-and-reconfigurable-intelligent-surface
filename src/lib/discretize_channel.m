function [discreteChannel] = discretize_channel(nInputs, nOutputs, symbolRatio, receivedPower, threshold)
	% Function:
	%	- obtain discrete channel based on channel probability distribution and bin boundaries
	%	- construct equivalent DMTC based on decision thresholds
    %
    % Input:
    %   - nInputs: size of input alphabet
	%	- nOutputs: number of output bins (or size of output alphabet)
	%	- symbolRatio: the duration ratio of the secondary symbol period over the primary symbol period
	%	- receivedPower: sorted received power per primary symbol under all tag input combinations
	%	- threshold: bin boundaries (or decision thresholds)
    %
    % Output:
	%	- discreteChannel: the DMC probability mass function after quantization (or the DMTC based on decision thresholds)
    %
    % Comment:
    %   - for a given tag input combination, the continous-output channel follows Erlang distribution
	%	- the discrete channel depends on bin boundaries (or decision thresholds)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	discreteChannel = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			discreteChannel(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
end
