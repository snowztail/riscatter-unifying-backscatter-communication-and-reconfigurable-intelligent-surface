function [capacity, inputDistribution] = blahut_arimoto(forwardTransition, tolerance)
	% Function:
    %   - compute the capacity and optimal input distribution of discrete memoryless channels
    %
    % Input:
    %   - forwardTransition [nInputs * nOutputs]: the channel transition probability matrix
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
    %   - capacity: maximum reliable rate for the input discrete memoryless channel
	%	- inputDistribution: optimal probability mass function of input distribution
    %
    % Comment:
    %   - backwardTransition: detecting transition probability matrix
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 07


	% * Get data
	nInputs = size(forwardTransition, 1);
	nOutputs = size(forwardTransition, 2);

	% * Initialize input distribution and backward transition matrix
	inputDistribution = ones(nInputs, 1) / nInputs;

	backwardTransition = zeros(nOutputs, nInputs);
	for iOutput = 1 : nOutputs
		for iInput = 1 : nInputs
			backwardTransition(iOutput, iInput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) / (transpose(inputDistribution) * forwardTransition(:, iOutput));
		end
	end

	% * Alternatively update input distribution and backward transition matrix
	capacity = 0;
	isConverged = false;
	while ~isConverged
		for iInput = 1 : nInputs
			inputDistribution(iInput) = prod(backwardTransition(:, iInput) .^ transpose(forwardTransition(iInput, :)));
		end
		inputDistribution = inputDistribution / sum(inputDistribution);

		for iOutput = 1 : nOutputs
			for iInput = 1 : nInputs
				backwardTransition(iOutput, iInput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) / (transpose(inputDistribution) * forwardTransition(:, iOutput));
			end
		end

		% * Update mutual information
		mutualInformation = zeros(nInputs, nOutputs);
		for iInput = 1 : nInputs
			for iOutput = 1 : nOutputs
				mutualInformation(iInput, iOutput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) * log2(backwardTransition(iOutput, iInput) / inputDistribution(iInput));
			end
		end
		mutualInformation = sum(mutualInformation, 'all');

		% * Test convergence
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end
end
