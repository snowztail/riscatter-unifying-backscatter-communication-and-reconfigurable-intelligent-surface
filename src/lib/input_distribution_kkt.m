function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr, tolerance)
	% Function:
	%	- obtain the tag input distribution that satisfies the KKT conditions of maximizing weighted sum primary-backscatter rate
    %
    % Input:
	%	- nTags: number of tags
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
    %	- iteratively update the input distribution by coordinate descent
	%	- returns KKT solution but without optimality guarantee
	%	- the marginal information function associated with codeword c_m_k of tag k satisfies:
	%		- I_k(c_m_k; Z) = C for P_k(c_m_k) > 0
	%		- I_k(c_m_k; Z) <= C for P_k(c_m_k) = 0
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Declare default tolerance
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
		tolerance = eps;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Initialization
	inputDistribution = ones(nTags, nStates) ./ nStates;
	combinationDistribution = combination_distribution(inputDistribution);
	equivalentDistribution = prod(combinationDistribution, 1);
	informationFunction = information_function(weight, symbolRatio, snr, equivalentDistribution, dmtc);
	marginalInformation = marginal_information(combinationDistribution, informationFunction);
	mutualInformation = equivalentDistribution * informationFunction;

	% * Iteratively update input distribution for all tags
	weightedSumRate = mutualInformation;
	isConverged = false;
	while ~isConverged
		% * Update input distribution, information functions associated with each codeword, marginal information of each codeword, and mutual information for each tag
		for iTag = 1 : nTags
			inputDistribution(iTag, :) = input_distribution_local(inputDistribution(iTag, :), marginalInformation(:, iTag));
			combinationDistribution = combination_distribution(inputDistribution);
			equivalentDistribution = prod(combinationDistribution, 1);
			informationFunction = information_function(weight, symbolRatio, snr, equivalentDistribution, dmtc);
			marginalInformation = marginal_information(combinationDistribution, informationFunction);
			mutualInformation = equivalentDistribution * informationFunction;
		end
		isConverged = abs(mutualInformation - weightedSumRate) <= tolerance;
		weightedSumRate = mutualInformation;
	end
end


function [inputDistribution] = input_distribution_local(inputDistribution, marginalInformation)
	nStates = size(inputDistribution, 2);
	inputDistribution_ = inputDistribution;
	for iState = 1 : nStates
		inputDistribution_(iState) = inputDistribution(iState) * exp(marginalInformation(iState)) / (inputDistribution * exp(marginalInformation));
	end
	inputDistribution = inputDistribution_;
end

function [marginalInformation] = marginal_information(combinationDistribution, informationFunction)
	[nTags, nInputs] = size(combinationDistribution);
	nStates = nthroot(nInputs, nTags);
	tagSet = transpose(1 : nTags);
	indexCombination = index_combination(nTags, nStates);
	marginalInformation = zeros(nTags, nStates, nInputs);
	for iTag = 1 : nTags
		for iState = 1 : nStates
			marginalSet = find(indexCombination(iTag, :) == iState);
			for iInput = marginalSet
				marginalInformation(iTag, iState, iInput) = prod(combinationDistribution(setdiff(tagSet, iTag), iInput), 1) * informationFunction(iInput);
			end
		end
	end
	marginalInformation = transpose(sum(marginalInformation, 3));
end
