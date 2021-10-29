function [capacity, inputDistribution] = blahut_arimoto_iterative(forwardTransition, tolerance)
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
    %   - iteratively update the input distribution and capacity
	%	- I(x; Y) = C for P(x) > 0
	%	- I(x; Y) <= C for P(x) = 0
    %
    % Author & Date: Yang (i@snowztail.com), 21 Oct 29

	% * Set default tolerance
	arguments
		forwardTransition;
		tolerance = eps;
	end

	% * Get data
	nInputs = size(forwardTransition, 1);
	nOutputs = size(forwardTransition, 2);

	% * Initialize input distribution and mutual information associated with each codeword
	inputDistribution = ones(nInputs, 1) / nInputs;
	mutualInformation = mutual_information(nInputs, nOutputs, inputDistribution, forwardTransition);

	% * Iteratively update the input distribution
	capacity = 0;
	isConverged = false;
	while ~isConverged
		inputDistribution = input_distribution(nInputs, inputDistribution, mutualInformation);

		% * Update mutual information
		mutualInformation = mutual_information(nInputs, nOutputs, inputDistribution, forwardTransition);

		% * Test convergence
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end
end


function [inputDistribution] = input_distribution(nInputs, inputDistribution, mutualInformation)
	for iInput = 1 : nInputs
		inputDistribution(iInput) = inputDistribution(iInput) * exp(mutualInformation(iInput)) / (transpose(inputDistribution) * exp(mutualInformation));
	end
end

function [mutualInformation] = mutual_information(nInputs, nOutputs, inputDistribution, forwardTransition)
	mutualInformation = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			mutualInformation(iInput, iOutput) = forwardTransition(iInput, iOutput) * log2(forwardTransition(iInput, iOutput) / (transpose(inputDistribution) * forwardTransition(:, iOutput)));
		end
	end
	mutualInformation = sum(mutualInformation, 2);
end
