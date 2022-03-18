function [rate, inputDistribution, equivalentDistribution] = alternating_optimization(weight, nTags, symbolRatio, snr, thresholdCandidate, dmc, receivedPower, equivalentDistribution, tolerance, options)
	% Function:
    %	- maximize the weighted sum of primary and total backscatter rate by alternatively optimizing input distribution and energy detection threshold
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- nTags: number of tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- thresholdCandidate [1 * (nBins + 1)]: candidate threshold values
    %	- dmc [(nStates ^ nTags) * nBins]: the transition probability matrix of the backscatter discrete memoryless MAC obtained by quantization
	%	- receivedPower [(nStates ^ nTags) * 1]: received power per primary symbol corresponding to each input letter combination combination
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- rate [1 * 2]: the achievable primary rate (nats per second per Hertz) and backscatter sum rate (nats per channel use)
    %
    % Comment:
	%	- alternatively update input (choose from 'cooperation', 'exhaustion', 'kkt', 'sca', 'marginalization', 'decomposition', 'randomization') and threshold (choose from 'smawk', 'dp', 'bisection', 'ml')
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 18

	% * Declare options for input and threshold design
	arguments
		weight;
		nTags;
		symbolRatio;
		snr;
		thresholdCandidate;
		dmc;
		receivedPower;
		equivalentDistribution = ones(1, size(snr, 1)) / size(snr, 1);
		tolerance = 1e-6;
		options.Input {mustBeMember(options.Input, ['cooperation', 'exhaustion', 'kkt', 'sca', 'marginalization', 'decomposition', 'randomization'])};
		options.Threshold {mustBeMember(options.Threshold, ['smawk', 'dp', 'bisection', 'ml'])};
	end

	% * Alternatively update threshold and input distribution until convergence
	weightedSumRate = 0;
	isConverged = false;
	while ~isConverged
		weightedSumRate_ = weightedSumRate;
		% * Threshold design
		switch options.Threshold
		case 'smawk'
			[~, dmtc] = threshold_smawk(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
		case 'dp'
			[~, dmtc] = threshold_dp(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
		case 'bisection'
			[~, dmtc] = threshold_bisection(thresholdCandidate, dmc, equivalentDistribution, receivedPower, symbolRatio);
		case 'ml'
			[~, dmtc] = threshold_ml(equivalentDistribution, receivedPower, symbolRatio);
		end
		% * Input design
		switch options.Input
		case 'cooperation'
			[inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
		case 'exhaustion'
			[inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_exhaustion(nTags, dmtc, weight, symbolRatio, snr);
		case 'kkt'
			[inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_kkt(nTags, dmtc, weight, symbolRatio, snr);
		case 'sca'
			[inputDistribution, equivalentDistribution, weightedSumRate] = input_distribution_sca(nTags, dmtc, weight, symbolRatio, snr);
		case 'marginalization'
			jointDistribution = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
			[inputDistribution, equivalentDistribution, weightedSumRate] = recovery_marginalization(jointDistribution, dmtc, weight, symbolRatio, snr);
		case 'decomposition'
			jointDistribution = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
			[inputDistribution, equivalentDistribution, weightedSumRate] = recovery_decomposition(jointDistribution, dmtc, weight, symbolRatio, snr);
		case 'randomization'
			jointDistribution = input_distribution_cooperation(nTags, dmtc, weight, symbolRatio, snr);
			[inputDistribution, equivalentDistribution, weightedSumRate] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr);
		end
		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
	[~, rate(1), rate(2)] = rate_weighted_sum(weight, symbolRatio, snr, equivalentDistribution, dmtc);
end
