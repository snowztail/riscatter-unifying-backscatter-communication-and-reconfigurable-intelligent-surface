function [rate, inputDistribution, threshold, beamformer] = block_coordinate_descent(weight, nTags, symbolRatio, equivalentChannel, txPower, noisePower, nBins, tolerance, options)
	% Function:
    %	- maximize the weighted sum of primary and total backscatter rate by alternatively optimizing input distribution and energy detection threshold
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- nTags: number of tags
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- txPower: average transmit power at the AP
	%	- noisePower: average noise power at the user
	%	- nBins: number of discretization bins over received signal
    %   - tolerance: minimum rate gain per iteration
	%	- options:
	%		- Input: choose from 'cooperation', 'exhaustion', 'kkt', 'sca'
	%		- Threshold: choose from 'smawk', 'dp', 'bisection', 'ml'
	%		- Recovery: choose from 'marginalization', 'decomposition', 'randomization'
    %
    % Output:
	%	- rate [2 * 1]:
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- threshold [1 * (nOutputs + 1)]: boundaries of decision regions
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %
    % Comment:
	%	- first perform non-coherent energy detection for backscatter links, then model their contribution within equivalent channel and decode primary message
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 18

	% * Declare options for input and threshold design
	arguments
		weight;
		nTags;
		symbolRatio;
		equivalentChannel;
		txPower;
		noisePower;
		nBins;
		tolerance = 1e-6;
		options.Input {mustBeMember(options.Input, ['cooperation', 'exhaustion', 'kkt', 'sca'])};
		options.Threshold {mustBeMember(options.Threshold, ['smawk', 'dp', 'bisection', 'ml'])};
		options.Recovery {mustBeMember(options.Recovery, ['marginalization', 'decomposition', 'randomization'])};
	end

	% * Get data
	nStates = nthroot(size(equivalentChannel, 1), nTags);

	% * Initialize input distribution, transmit beamformer, and decision threshold
	inputDistribution = ones(nTags, nStates) / nStates;
	equivalentDistribution = prod(combination_distribution(inputDistribution), 1);
	beamformer = sqrt(txPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
	[threshold, dmtc] = threshold_ml(symbolRatio, equivalentChannel, noisePower, beamformer);

	% * Iteratively update input distribution, tranmit beamformer, and decision threshold
	weightedSumRate = 0;
	isConverged = false;
	while ~isConverged
		weightedSumRate_ = weightedSumRate;

		% * Input distribution design
		switch options.Input
		case 'exhaustion'
			[inputDistribution, equivalentDistribution] = input_distribution_exhaustion(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc);
		case 'kkt'
			[inputDistribution, equivalentDistribution] = input_distribution_kkt(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc);
		case 'sca'
			[inputDistribution, equivalentDistribution] = input_distribution_sca(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc, tolerance);
		case 'cooperation'
			[jointDistribution, equivalentDistribution] = input_distribution_cooperation(weight, nTags, symbolRatio, equivalentChannel, noisePower, beamformer, dmtc);
			if isfield(options, 'Recovery')
				switch options.Recovery
				case 'marginalization'
					[inputDistribution, equivalentDistribution] = recovery_marginalization(weight, symbolRatio, equivalentChannel, noisePower, jointDistribution, beamformer, dmtc);
				case 'decomposition'
					[inputDistribution, equivalentDistribution] = recovery_decomposition(weight, symbolRatio, equivalentChannel, noisePower, jointDistribution, beamformer, dmtc);
				case 'randomization'
					[inputDistribution, equivalentDistribution] = recovery_randomization(weight, symbolRatio, equivalentChannel, noisePower, jointDistribution, beamformer, dmtc);
				end
			end
		end

		% * Beamforming design
		beamformer = beamformer_sca(weight, symbolRatio, equivalentChannel, txPower, noisePower, equivalentDistribution, threshold, tolerance);

		% * Threshold design
		switch options.Threshold
		case 'smawk'
			[threshold, dmtc] = threshold_smawk(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer);
		case 'dp'
			[threshold, dmtc] = threshold_dp(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer);
		case 'bisection'
			[threshold, dmtc] = threshold_bisection(symbolRatio, equivalentChannel, noisePower, nBins, equivalentDistribution, beamformer);
		case 'ml'
			[threshold, dmtc] = threshold_ml(symbolRatio, equivalentChannel, noisePower, beamformer);
		end

		% * Test convergence
		[weightedSumRate, rate] = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end
end
