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
	%	- inputDistribution: optimal probability mass function of input distribution (Lagrange multiplier)
    %
    % Comment:
    %   - backwardTransition: detecting transition probability matrix (Bayes' Theorem)
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 07

	% * Set default tolerance
	arguments
		forwardTransition;
		tolerance = eps;
	end

	% * Get data
	nInputs = size(forwardTransition, 1);
	nOutputs = size(forwardTransition, 2);

	% * Initialize input distribution and backward transition matrix
	inputDistribution = ones(nInputs, 1) / nInputs;
	backwardTransition = backward_transition(nInputs, nOutputs, inputDistribution, forwardTransition);

	% * Alternatively update input distribution and backward transition matrix
	capacity = 0;
	isConverged = false;
	while ~isConverged
		inputDistribution = input_distribution(nInputs, forwardTransition, backwardTransition);
		backwardTransition = backward_transition(nInputs, nOutputs, inputDistribution, forwardTransition);

		% * Update mutual information
		mutualInformation = mutual_information(nInputs, nOutputs, forwardTransition, backwardTransition, inputDistribution);

		% * Test convergence
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end
end


function [backwardTransition] = backward_transition(nInputs, nOutputs, inputDistribution, forwardTransition)
	backwardTransition = zeros(nOutputs, nInputs);
	for iOutput = 1 : nOutputs
		for iInput = 1 : nInputs
			backwardTransition(iOutput, iInput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) / (transpose(inputDistribution) * forwardTransition(:, iOutput));
		end
	end
end

function [inputDistribution] = input_distribution(nInputs, forwardTransition, backwardTransition)
	inputDistribution = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		inputDistribution(iInput) = prod(backwardTransition(:, iInput) .^ transpose(forwardTransition(iInput, :)));
	end
	inputDistribution = inputDistribution / sum(inputDistribution);
end

function [mutualInformation] = mutual_information(nInputs, nOutputs, forwardTransition, backwardTransition, inputDistribution)
	mutualInformation = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			mutualInformation(iInput, iOutput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) * log2(backwardTransition(iOutput, iInput) / inputDistribution(iInput));
		end
	end
	mutualInformation = sum(mutualInformation, 'all');
end
