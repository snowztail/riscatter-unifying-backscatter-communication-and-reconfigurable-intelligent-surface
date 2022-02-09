function [wsr, primaryRate, secondaryRate, inputDistribution, equivalentDistribution] = input_distribution(weight, symbolRatio, snr, dmtc, nTags, tolerance)
	% Function:
    %   - compute the weighted sum rate of the primary user and all backscatter tags
	%	- obtain the optimal tag input distribution for a given discrete memoryless MAC
    %
    % Input:
	%	- weight: the priority of the primary link
	%	- symbolRatio: the duration ratio of the secondary symbol period over the primary symbol period
	%	- snr [nInputs * 1]: signal-to-noise ratio of the primary link corresponding to all tag input combinations
    %   - dmtc [nInputs * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- nTags: number of tags
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- wsr: weighted sum rate of the primary user and all backscatter tags
	%	- primaryRate: the achievable rate for the primary link (bps/Hz)
	%	- secondaryRate: the achievable rate for the secondary link (bpcu)
	%	- inputDistribution: optimal input probability distribution
	%	- equivalentDistribution: optimal input combination probability distribution
    %
    % Comment:
    %   - iteratively update the input distribution and achievable rates
	%	- the optimal marginal information function associated with codeword c_m_k of tag k satisfies:
	%		- I_k(c_m_k; Z) = C for P_k(c_m_k) > 0
	%		- I_k(c_m_k; Z) <= C for P_k(c_m_k) = 0
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Set default tolerance
	arguments
		weight;
		symbolRatio;
		snr;
		dmtc;
		nTags;
		tolerance = 1e-12;
	end

	% * Ensure non-zero transitional probability as required by Blahut-Arimoto algorithm
	dmtc(dmtc < eps) = eps;

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Initialize
	inputDistribution = normr(sqrt(rand(nTags, nStates))) .^ 2;
	combinationDistribution = combination_distribution(inputDistribution);
	equivalentDistribution = prod(combinationDistribution, 1);
	primaryInformationFunction = primary_information_function(symbolRatio, snr);
	secondaryInformationFunction = secondary_information_function(equivalentDistribution, dmtc);
	primaryRate = equivalentDistribution * primaryInformationFunction;
	secondaryRate = equivalentDistribution * secondaryInformationFunction;
	informationFunction = weight * primaryInformationFunction + (1 - weight) * secondaryInformationFunction;
	marginalInformation = marginal_information(combinationDistribution, informationFunction);
	mutualInformation = equivalentDistribution * informationFunction;

	% * Iteratively update input distribution for all tags
	wsr = mutualInformation;
	isConverged = false;
	while ~isConverged
		% * Update input distribution, information functions associated with each codeword, marginal information of each codeword, and mutual information for each tag
		for iTag = 1 : nTags
			inputDistribution(iTag, :) = input_distribution_local(inputDistribution(iTag, :), marginalInformation(:, iTag));
			combinationDistribution = combination_distribution(inputDistribution);
			equivalentDistribution = prod(combinationDistribution, 1);
			primaryInformationFunction = primary_information_function(symbolRatio, snr);
			secondaryInformationFunction = secondary_information_function(equivalentDistribution, dmtc);
			primaryRate = equivalentDistribution * primaryInformationFunction;
			secondaryRate = equivalentDistribution * secondaryInformationFunction;
			informationFunction = weight * primaryInformationFunction + (1 - weight) * secondaryInformationFunction;
			marginalInformation = marginal_information(combinationDistribution, informationFunction);
			mutualInformation = equivalentDistribution * informationFunction;
		end
		isConverged = abs(mutualInformation - wsr) <= tolerance;
		wsr = mutualInformation;
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

function [combinationDistribution] = combination_distribution(inputDistribution)
	[nTags, nStates] = size(inputDistribution);
	nInputs = nStates ^ nTags;
	tagSet = transpose(1 : nTags);
	indexCombination = nested_combvec(1 : nStates, nTags);
	combinationDistribution = zeros(nTags, nInputs);
	for iInput = 1 : nInputs
		combinationDistribution(:, iInput) = inputDistribution(sub2ind(size(inputDistribution), tagSet, indexCombination(:, iInput)));
	end
end

function [primaryInformationFunction] = primary_information_function(symbolRatio, snr)
	nInputs = size(snr, 1);
	primaryInformationFunction = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		primaryInformationFunction(iInput) = symbolRatio * log2(1 + snr(iInput));
	end
end

function [secondaryInformationFunction] = secondary_information_function(equivalentDistribution, dmtc)
	[nInputs, nOutputs] = size(dmtc);
	secondaryInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			secondaryInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log2(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	secondaryInformationFunction = sum(secondaryInformationFunction, 2);
end

function [marginalInformation] = marginal_information(combinationDistribution, informationFunction)
	[nTags, nInputs] = size(combinationDistribution);
	nStates = nthroot(nInputs, nTags);
	tagSet = transpose(1 : nTags);
	indexCombination = nested_combvec(1 : nStates, nTags);
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
