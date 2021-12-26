function [capacity, inputDistribution] = multiuser_blahut_arimoto(dmc, nUsers, tolerance)
	% Function:
    %   - compute the sum capacity and optimal multiuser input distribution of discrete memoryless MAC
    %
    % Input:
    %   - dmc [nInputs * nOutputs]: the transition probability matrix of discrete memoryless MAC
	%	- nUsers: number of users
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
    %   - capacity: maximum achievable sum rate of the discrete memoryless MAC
	%	- inputDistribution: optimal multiuser input probability distribution corresponding to the alphabets
    %
    % Comment:
    %   - iteratively update the input distribution and capacity
	%	- the optimal marginal information function associated with codeword x of user k satisfies:
	%		- I_k(x; Y) = C for P_k(x) > 0
	%		- I_k(x; Y) <= C for P_k(x) = 0
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 21 Dec 26

	% * Set default tolerance
	arguments
		dmc;
		nUsers;
		tolerance = eps;
	end

	% * Ensure non-zero transitional probability as required
	dmc(dmc < eps) = eps;

	% * Get data
	nInputs = size(dmc, 1);
	nOutputs = size(dmc, 2);
	nStates = nthroot(nInputs, nUsers);

	% * Initialize
	inputDistribution = normr(sqrt(rand(nUsers, nStates))) .^ 2;
	combinationDistribution = combination_distribution(nUsers, nStates, inputDistribution);
	informationFunction = information_function(nInputs, nOutputs, prod(combinationDistribution), dmc);
	marginalInformation = marginal_information(nUsers, nStates, combinationDistribution, informationFunction);
	mutualInformation = prod(combinationDistribution) * informationFunction;

	% * Iteratively update input distribution, information function associated with each codeword, marginal information of each codeword of each user, and mutual information
	capacity = mutualInformation;
	isConverged = false;
	while ~isConverged
		for iUser = 1 : nUsers
			inputDistribution_ = zeros(1, nStates);
			for iState = 1 : nStates
				inputDistribution_(iState) = inputDistribution(iUser, iState) * exp(marginalInformation(iState, iUser)) / (inputDistribution(iUser, :) * exp(marginalInformation(:, iUser)));
			end
			inputDistribution(iUser, :) = inputDistribution_;
			combinationDistribution = combination_distribution(nUsers, nStates, inputDistribution);
			informationFunction = information_function(nInputs, nOutputs, prod(combinationDistribution), dmc);
			marginalInformation = marginal_information(nUsers, nStates, combinationDistribution, informationFunction);
			mutualInformation = prod(combinationDistribution) * informationFunction;
		end
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end

	% * Discard codewords with negligible probability
	inputDistribution(inputDistribution < eps) = 0;
end


function [inputDistribution] = input_distribution(nInputs, inputDistribution, informationFunction)
	for iInput = 1 : nInputs
		inputDistribution(iInput) = inputDistribution(iInput) * exp(informationFunction(iInput)) / (inputDistribution * exp(informationFunction));
	end
end

function [combinationDistribution] = combination_distribution(nUsers, nStates, inputDistribution)
	nInputs = nStates ^ nUsers;
	userSet = transpose(1 : nUsers);
	indexCombination = nested_combvec(1 : nStates, nUsers);
	combinationDistribution = zeros(nUsers, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), userSet, indexCombination(:, iInput)));
	end
end

function [informationFunction] = information_function(nInputs, nOutputs, inputDistribution, dmc)
	informationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			informationFunction(iInput, iOutput) = dmc(iInput, iOutput) * log2(dmc(iInput, iOutput) / (inputDistribution * dmc(:, iOutput)));
		end
	end
	informationFunction = sum(informationFunction, 2);
end

function [marginalInformation] = marginal_information(nUsers, nStates, combinationDistribution, informationFunction)
	nMarginals = nStates ^ (nUsers - 1);
	userSet = transpose(1 : nUsers);
	marginalInformation = zeros(nUsers, nStates, nMarginals);
	for iUser = 1 : nUsers
		for iState = 1 : nStates
			for iMarginal = 1 : nMarginals
				iInput = (iState - 1) * nMarginals + iMarginal;
				marginalInformation(iUser, iState, iMarginal) = prod(combinationDistribution(setdiff(userSet, iUser), iInput)) * informationFunction(iInput);
			end
		end
	end
	marginalInformation = transpose(sum(marginalInformation, 3));
end
