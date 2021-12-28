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
	nStates = nthroot(size(dmc, 1), nUsers);

	% * Initialize
	inputDistribution = normr(sqrt(rand(nUsers, nStates))) .^ 2;
	combinationDistribution = combination_distribution(inputDistribution);
	jointDistribution = prod(combinationDistribution, 1);
	informationFunction = information_function(jointDistribution, dmc);
	marginalInformation = marginal_information(combinationDistribution, informationFunction);
	mutualInformation = jointDistribution * informationFunction;

	% * Iteratively update input distribution for all users
	capacity = mutualInformation;
	isConverged = false;
	while ~isConverged
		% * Update input distribution, information function associated with each codeword, marginal information of each codeword, and mutual information for a single user
		for iUser = 1 : nUsers
			inputDistribution(iUser, :) = input_distribution(inputDistribution(iUser, :), marginalInformation(:, iUser));
			combinationDistribution = combination_distribution(inputDistribution);
			jointDistribution = prod(combinationDistribution, 1);
			informationFunction = information_function(jointDistribution, dmc);
			marginalInformation = marginal_information(combinationDistribution, informationFunction);
			mutualInformation = jointDistribution * informationFunction;
		end
		isConverged = abs(mutualInformation - capacity) <= tolerance;
		capacity = mutualInformation;
	end

	% * Discard codewords with negligible probability
	inputDistribution(inputDistribution < eps) = 0;
end


function [inputDistribution] = input_distribution(inputDistribution, marginalInformation)
	nStates = size(inputDistribution, 2);
	inputDistribution_ = inputDistribution;
	for iState = 1 : nStates
		inputDistribution_(iState) = inputDistribution(iState) * exp(marginalInformation(iState)) / (inputDistribution * exp(marginalInformation));
	end
	inputDistribution = inputDistribution_;
end

function [combinationDistribution] = combination_distribution(inputDistribution)
	[nUsers, nStates] = size(inputDistribution);
	nInputs = nStates ^ nUsers;
	userSet = transpose(1 : nUsers);
	indexCombination = nested_combvec(1 : nStates, nUsers);
	combinationDistribution = zeros(nUsers, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), userSet, indexCombination(:, iInput)));
	end
end

function [informationFunction] = information_function(jointDistribution, dmc)
	[nInputs, nOutputs] = size(dmc);
	informationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			informationFunction(iInput, iOutput) = dmc(iInput, iOutput) * log2(dmc(iInput, iOutput) / (jointDistribution * dmc(:, iOutput)));
		end
	end
	informationFunction = sum(informationFunction, 2);
end

function [marginalInformation] = marginal_information(combinationDistribution, informationFunction)
	[nUsers, nInputs] = size(combinationDistribution);
	nStates = nthroot(nInputs, nUsers);
	userSet = transpose(1 : nUsers);
	indexCombination = nested_combvec(1 : nStates, nUsers);
	marginalInformation = zeros(nUsers, nStates, nInputs);
	for iUser = 1 : nUsers
		for iState = 1 : nStates
			marginalSet = find(indexCombination(iUser, :) == iState);
			% marginalSet = transpose(vec((0 : nStates ^ (iUser - 1) - 1) * nStates ^ (nUsers - iUser + 1) + transpose((iState - 1) * nStates ^ (nUsers - iUser) + 1 : iState * nStates ^ (nUsers - iUser))));
			for iInput = marginalSet
				marginalInformation(iUser, iState, iInput) = prod(combinationDistribution(setdiff(userSet, iUser), iInput), 1) * informationFunction(iInput);
			end
		end
	end
	marginalInformation = transpose(sum(marginalInformation, 3));
end
