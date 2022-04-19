function [inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_kkt(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc, tolerance)
	% Function:
	%	- obtain the tag input distribution that satisfies the KKT conditions of maximizing weighted sum primary-backscatter rate
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- nTags: number of tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
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
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_K)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Declare default tolerance
	arguments
		weight;
		nTags;
		symbolRatio;
		equivalentChannel;
		noisePower;
		beamformer;
		dmtc;
		tolerance = 1e-6;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Initialize input distribution as uniform distribution
	inputDistribution = ones(nTags, nStates) / nStates;
	combinationDistribution = combination_distribution(inputDistribution);
	equivalentDistribution = prod(combinationDistribution, 1);
	informationFunction = information_function(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
	marginalInformation = marginal_information(combinationDistribution, informationFunction);
	weightedSumRate = equivalentDistribution * informationFunction;

	% * Iteratively update input distribution for all tags
	isConverged = false;
	while ~isConverged
		weightedSumRate_ = weightedSumRate;
		% * Update input distribution, information functions associated with each codeword, marginal information of each codeword, and mutual information for each tag
		for iTag = 1 : nTags
			inputDistribution(iTag, :) = input_distribution(inputDistribution(iTag, :), marginalInformation(:, iTag));
			combinationDistribution = combination_distribution(inputDistribution);
			equivalentDistribution = prod(combinationDistribution, 1);
			informationFunction = information_function(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
			marginalInformation = marginal_information(combinationDistribution, informationFunction);
			weightedSumRate = equivalentDistribution * informationFunction;
		end
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
end


function [inputDistribution] = input_distribution(inputDistribution, marginalInformation)
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
