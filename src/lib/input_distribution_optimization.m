function [inputDistribution, equivalentDistribution, weightedSumRateUpperBound, weightedSumRateLowerBound] = input_distribution_optimization(nTags, dmtc, weight, symbolRatio, snr)
	% Function:
	%	- optimize the joint input distribution of all tags (corresponding to the optimal input with full transmit cooperation)
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by randomization
    %
    % Input:
	%	- nTags: number of tags
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRateUpperBound: maximum achievable weighted sum rate with tag transmit correlation
	%	- weightedSumRateLowerBound: maximum achievable weighted sum rate without tag transmit correlation
    %
    % Comment:
	%	- joint input optimization
	%		- tags are assumed with equal weight (hence sum rate of tags is considered)
	%	- rank-1 approximation
    %		- generate random samples of candidate tag probability vectors whose outer product expectation equals to the joint input distribution array
	%		- equivalent to a random search with guiduance on the correlation matrix of the distributions
	%		- the number of kernels and random samples are designable (performance-complexity tradeoff)
	%		- the candidates follows uniform distribution within a sphere bounded within the probability simplex
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 23

	% * Declare default tolerance
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Optimize joint input probability distribution of all tags
	cvx_begin
		variables jointDistribution(nStates * ones(1, nTags));
		equivalentDistribution = transpose(vec(permute(jointDistribution, nTags : -1 : 1)));
		primaryRate = equivalentDistribution * information_function_primary(symbolRatio, snr);
		backscatterRate = backscatter_rate(equivalentDistribution, dmtc);
		weightedSumRateUpperBound = [primaryRate, backscatterRate] * weight;

		maximize weightedSumRateUpperBound
		subject to
			jointDistribution == nonnegative(size(jointDistribution));
			sum(equivalentDistribution) == 1;
	cvx_end
	[inputDistribution, equivalentDistribution, weightedSumRateLowerBound] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr);
end


function [backscatterRate] = backscatter_rate(equivalentDistribution, dmtc)
	nOutputs = size(dmtc, 2);
	backscatterRate = cvx(zeros(nOutputs, 1));
	for iOutput = 1 : nOutputs
		backscatterRate(iOutput) = entr(equivalentDistribution * dmtc(:, iOutput)) - equivalentDistribution * entr(dmtc(:, iOutput));
	end
	backscatterRate = sum(backscatterRate);
end
