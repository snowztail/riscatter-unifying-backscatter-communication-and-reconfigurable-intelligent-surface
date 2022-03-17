function [jointDistribution, equivalentDistribution, weightedSumRate] = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr, tolerance)
	% Function:
	%	- obtain the optimal joint input distribution with full transmit cooperation
    %
    % Input:
	%	- nTags: number of tags
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: maximum achievable weighted sum of primary rate and total backscatter rate with full tag transmit cooperation
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
	%	- with full transmit cooperation, the tags are equivalent to a single information source with input distribution in (nStates ^ nTags)-dimensional probability simplex
    %   - in general, tags encode independently and it is hard to approximate the joint distribution array without cooperation
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Declare default tolerance
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
		tolerance = eps;
	end

	% * Ensure non-zero transitional probability
	dmtc(dmtc < eps) = eps;

	% * Get data
	nInputs = size(dmtc, 1);
	nStates = nthroot(nInputs, nTags);

	% * Initialization
	equivalentDistribution = ones(1, nInputs) ./ nInputs;
	informationFunction = information_function(weight, symbolRatio, snr, equivalentDistribution, dmtc);
	weightedSumRate = equivalentDistribution * informationFunction;

	% * Iteratively update equivalent distribution, information function associated with each (joint) codeword, and weighted sum achievable rate
	isConverged = false;
	while ~isConverged
		weightedSumRate_ = weightedSumRate;
		equivalentDistribution = equivalent_distribution(equivalentDistribution, informationFunction);
		jointDistribution = permute(reshape(equivalentDistribution, nStates * ones(1, nTags)), nTags : - 1 : 1);
		informationFunction = information_function(weight, symbolRatio, snr, equivalentDistribution, dmtc);
		weightedSumRate = equivalentDistribution * informationFunction;
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
end


function [equivalentDistribution] = equivalent_distribution(equivalentDistribution, informationFunction)
	nInputs = size(equivalentDistribution, 2);
	equivalentDistribution_ = equivalentDistribution;
	for iInput = 1 : nInputs
		equivalentDistribution_(iInput) = equivalentDistribution(iInput) * exp(informationFunction(iInput)) / (equivalentDistribution * exp(informationFunction));
	end
	equivalentDistribution = equivalentDistribution_;
end
