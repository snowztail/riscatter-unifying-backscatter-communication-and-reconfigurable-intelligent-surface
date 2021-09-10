function [primaryRate, secondaryRate, inputDistribution] = weighted_sum_rate(forwardTransition, snr, sortedChannel, weight, tolerance)
	% Function:
    %   - compute the maximum weighted sum rate and optimal tag input distribution for symbiotic radio systems
    %
    % Input:
    %   - forwardTransition [nInputs * nOutputs]: the channel transition probability matrix
	%	- snr: average transmit SNR
	%	- sortedChannel [nStates * 1]: sorted equivalent channel in received energy ascending order
	%	- weight: the priority of the primary link
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- primaryRate: the achievable rate for the primary link (bps/Hz)
	%	- secondaryRate: the achievable rate for the secondary link (bpcu)
	%	- inputDistribution: optimal probability mass function of input distribution  (Lagrange multiplier)
    %
    % Comment:
    %   - backwardTransition: detecting transition probability matrix (Bayes' Theorem)
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 10


	% * Get data
	nInputs = size(forwardTransition, 1);
	nOutputs = size(forwardTransition, 2);

	% * Initialize input distribution and backward transition matrix
	inputDistribution = ones(nInputs, 1) / nInputs;
	backwardTransition = backward_transition(nInputs, nOutputs, inputDistribution, forwardTransition);

	% * Alternating optimization
	weightedSumCapacity = 0;
	isConverged = false;
	while ~isConverged
		% * Update input distribution and backward transition matrix
		inputDistribution = input_distribution(nInputs, forwardTransition, backwardTransition, snr, sortedChannel, weight);
		backwardTransition = backward_transition(nInputs, nOutputs, inputDistribution, forwardTransition);

		% * Update achievable rates
		primaryRate = primary_rate(nInputs, inputDistribution, snr, sortedChannel);
		secondaryRate = secondary_rate(nInputs, nOutputs, forwardTransition, backwardTransition, inputDistribution);
		weightedSumRate = weight * primaryRate + (1 - weight) * secondaryRate;

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumCapacity) <= tolerance || any(inputDistribution <= eps);
		weightedSumCapacity = weightedSumRate;
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

function [inputDistribution] = input_distribution(nInputs, forwardTransition, backwardTransition, snr, sortedChannel, weight)
	inputDistribution = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		inputDistribution(iInput) = (1 + snr * abs(sortedChannel(iInput)) ^ 2) ^ (weight / (1 - weight)) * prod(backwardTransition(:, iInput) .^ transpose(forwardTransition(iInput, :)));
	end
	inputDistribution = inputDistribution / sum(inputDistribution);
end

function [primaryRate] = primary_rate(nInputs, inputDistribution, snr, sortedChannel)
	primaryRate = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryRate(iInput) = inputDistribution(iInput) * log2(1 + snr * abs(sortedChannel(iInput)) ^ 2);
	end
	primaryRate = sum(primaryRate);
end

function [secondaryRate] = secondary_rate(nInputs, nOutputs, forwardTransition, backwardTransition, inputDistribution)
	secondaryRate = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			secondaryRate(iInput, iOutput) = inputDistribution(iInput) * forwardTransition(iInput, iOutput) * log2(backwardTransition(iOutput, iInput) / inputDistribution(iInput));
		end
	end
	secondaryRate = sum(secondaryRate, 'all');
end
