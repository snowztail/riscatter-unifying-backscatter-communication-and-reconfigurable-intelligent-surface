function [capacity, inputDistribution] = blahut_arimoto(dmc, tolerance)
	% Function:
    %   - compute the capacity and optimal input distribution of discrete memoryless channels
    %
    % Input:
    %   - dmc [nInputs * nOutputs]: the transition probability matrix of discrete memoryless channel
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
    %   - capacity: maximum achievable rate of DMC
	%	- inputDistribution: optimal input probability distribution corresponding to the alphabet
    %
    % Comment:
    %   - iteratively update the input distribution and capacity
	%	- the optimal information function associated with codeword x satisfies:
	%		- I(x; Y) = C for P(x) > 0
	%		- I(x; Y) <= C for P(x) = 0
    %
    % Author & Date: Yang (i@snowztail.com), 21 Oct 29

	% * Set default tolerance
	arguments
		dmc;
		tolerance = eps;
	end

	% * Ensure non-zero transitional probability as required
	dmc(dmc < eps) = eps;

	% * Get data
	nInputs = size(dmc, 1);

	% * Initialize
	inputDistribution = normr(sqrt(rand(1, nInputs))) .^ 2;
	informationFunction = information_function(inputDistribution, dmc);
	mutualInformation = inputDistribution * informationFunction;

	% * Iteratively update input distribution, information function associated with each codeword, and mutual information
	capacity = mutualInformation;
	isConverged = false;
	while ~isConverged
		inputDistribution = input_distribution(inputDistribution, informationFunction);
		informationFunction = information_function(inputDistribution, dmc);
		mutualInformation = inputDistribution * informationFunction;
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end

	% * Discard codewords with negligible probability
	inputDistribution(inputDistribution < eps) = 0;
end


function [inputDistribution] = input_distribution(inputDistribution, informationFunction)
	nInputs = size(inputDistribution, 2);
	inputDistribution_ = inputDistribution;
	for iInput = 1 : nInputs
		inputDistribution_(iInput) = inputDistribution(iInput) * exp(informationFunction(iInput)) / (inputDistribution * exp(informationFunction));
	end
	inputDistribution = inputDistribution_;
end

function [informationFunction] = information_function(inputDistribution, dmc)
	[nInputs, nOutputs] = size(dmc);
	informationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			informationFunction(iInput, iOutput) = dmc(iInput, iOutput) * log2(dmc(iInput, iOutput) / (inputDistribution * dmc(:, iOutput)));
		end
	end
	informationFunction = sum(informationFunction, 2);
end
