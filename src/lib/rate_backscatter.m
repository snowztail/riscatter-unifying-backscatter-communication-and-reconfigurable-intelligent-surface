function [backscatterRate] = rate_backscatter(equivalentDistribution, dmtc)
	% Function:
    %	- compute the total backscatter rate for a given input combination distribution and discrete thresholding MAC
    %
    % Input:
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
    %
    % Output:
	%	- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	% * Get data
	[nInputs, nOutputs] = size(dmtc);

	if isa(equivalentDistribution, 'cvx')
		% * Concave expression used by CVX
		backscatterRate = cvx(zeros(nOutputs, 1));
		for iOutput = 1 : nOutputs
			backscatterRate(iOutput) = entr(equivalentDistribution * dmtc(:, iOutput)) - equivalentDistribution * entr(dmtc(:, iOutput));
		end
		backscatterRate = sum(backscatterRate);
	else
		% * Numerical result
		backscatterInformationFunction = zeros(nInputs, nOutputs);
		for iInput = 1 : nInputs
			for iOutput = 1 : nOutputs
				backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
			end
		end
		backscatterRate = equivalentDistribution * sum(backscatterInformationFunction, 2);
	end
end
