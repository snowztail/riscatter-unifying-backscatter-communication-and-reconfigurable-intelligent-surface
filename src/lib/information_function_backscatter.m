function [backscatterInformationFunction] = information_function_backscatter(equivalentDistribution, dmtc)
	% Function:
    %	- compute the backscatter information function associated with each input letter combination for a given DMC and input distribution
    %
    % Input:
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- backscatterInformationFunction [(nStates ^ nTags) * 1]: backscatter information function associated with each input letter combination
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 19

	[nInputs, nOutputs] = size(dmtc);
	backscatterInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	backscatterInformationFunction = sum(backscatterInformationFunction, 2);
end
