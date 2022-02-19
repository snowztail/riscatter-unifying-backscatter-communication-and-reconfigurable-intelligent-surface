function [primaryInformationFunction] = information_function_primary(snr, symbolRatio)
	% Function:
    %   - compute the primary information function associated with each tag input letter combination
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to each input letter combination
    %
    % Output:
	%	- primaryInformationFunction [(nStates ^ nTags) * 1]: primary information function associated with each input letter
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 19

	nInputs = size(snr, 1);
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log2(1 + snr(iInput));
	end
end
