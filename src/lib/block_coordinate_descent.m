function [rate, distribution, threshold, beamforming] = block_coordinate_descent(nTags, symbolRatio, transmitPower, noisePower, nBins, weight, equivalentChannel, tolerance, Options)
	% Function:
	%	- iteratively update the input distribution, detection threshold, and beamforming vector to maximize weighted sum rate
    %
    % Input:
	%	- nTags: number of tags
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- nBins: number of receive energy quantization bins
	%	- weight: relative priority of primary link
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
    %   - tolerance: minimum rate gain ratio per iteration
	%	- Options: algorithm options as structure of name-value pairs
    %
    % Output:
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
	%	- beamforming [nTxs x 1]: transmit beamforming vector
    %
    % Comment:
	%	- user first jointly decode all tags by energy detection, then model backscatter contribution within equivalent channel
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 18

	% * Set default tolerance
	arguments
		nTags;
		symbolRatio;
		transmitPower;
		noisePower;
		nBins;
		weight;
		equivalentChannel;
		tolerance = 1e-6;
		Options.Distribution {mustBeMember(Options.Distribution, ['exhaustion', 'kkt', 'sca', 'cooperation'])};
		Options.Threshold {mustBeMember(Options.Threshold, ['smawk', 'dp', 'bisection', 'ml'])};
		Options.Beamforming {mustBeMember(Options.Beamforming, ['pgd', 'emrt', 'dmrt'])};
		Options.Recovery {mustBeMember(Options.Recovery, ['marginalization', 'decomposition', 'randomization'])};
	end

	% * Get data
	nInputs = size(equivalentChannel, 2);
	nStates = nthroot(nInputs, nTags);
	directChannel = sum(equivalentChannel, 2);

	% ! Initialize distribution and beamforming by previous solutions
	persistent Initializer

	% * No previous solutions, use uniform distribution and MRT towards direct/ergodic channels (equal under uniform distribution)
	if isempty(Initializer)
		Initializer.distribution = normalize(ones(nStates, nTags), 'norm', 1);
		Initializer.beamforming = sqrt(transmitPower) * directChannel / norm(directChannel);
	end

	% * Apply initializer
	distribution = Initializer.distribution;
	beamforming = Initializer.beamforming;
	equivalentDistribution = prod(tuple_tag(distribution), 2);

	% * Initialize decision threshold
	receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
	snr = receivePower / noisePower;
	thresholdDomain = domain_threshold(symbolRatio, nBins, receivePower);
	binDmc = dmc_integration(symbolRatio, receivePower, thresholdDomain);
	switch Options.Threshold
	case 'smawk'
		threshold = threshold_smawk(equivalentDistribution, thresholdDomain, binDmc);
	case 'dp'
		threshold = threshold_dp(equivalentDistribution, thresholdDomain, binDmc);
	case 'bisection'
		threshold = threshold_bisection(symbolRatio, receivePower, equivalentDistribution, thresholdDomain, binDmc);
	case 'ml'
		threshold = threshold_ml(symbolRatio, receivePower);
	end

	% * Construct DMTC and recover i/o mapping
	dmac = dmc_integration(symbolRatio, receivePower, threshold);
	[~, sortIndex] = sort(receivePower);
	dmac(:, sortIndex) = dmac;

	% * Block coordinate descent
	wsr = rate_weighted(weight, snr, equivalentDistribution, dmac);
	isConverged = false;

	while ~isConverged
		% * Update iteration index
		wsr_ = wsr;

		% * Input probability distribution
		switch Options.Distribution
		case 'exhaustion'
			[distribution, equivalentDistribution] = distribution_exhaustion(nTags, weight, snr, dmac);
		case 'kkt'
			[distribution, equivalentDistribution] = distribution_kkt(nTags, weight, snr, dmac);
		case 'sca'
			[distribution, equivalentDistribution] = distribution_sca(nTags, weight, snr, dmac);
		case 'cooperation'
			[jointDistribution, equivalentDistribution] = distribution_cooperation(nTags, weight, snr, dmac);
			if ~isempty(Options.Recovery)
				switch Options.Recovery
				case 'marginalization'
					[distribution, equivalentDistribution] = recovery_marginalization(jointDistribution);
				case 'decomposition'
					[distribution, equivalentDistribution] = recovery_decomposition(jointDistribution);
				case 'randomization'
					[distribution, equivalentDistribution] = recovery_randomization(weight, snr, jointDistribution, dmac);
				end
			end
		end

		% * Transmit beamforming
		switch Options.Beamforming
		case 'pgd'
			beamforming = beamforming_pgd(symbolRatio, weight, transmitPower, noisePower, equivalentChannel, equivalentDistribution, threshold);
		case 'emrt'
			beamforming = beamforming_emrt(transmitPower, equivalentChannel, equivalentDistribution);
		case 'dmrt'
			beamforming = beamforming_dmrt(transmitPower, directChannel);
		end

		% * Decision threshold
		receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
		snr = receivePower / noisePower;
		thresholdDomain = domain_threshold(symbolRatio, nBins, receivePower);
		binDmc = dmc_integration(symbolRatio, receivePower, thresholdDomain);
		switch Options.Threshold
		case 'smawk'
			threshold = threshold_smawk(equivalentDistribution, thresholdDomain, binDmc);
		case 'dp'
			threshold = threshold_dp(equivalentDistribution, thresholdDomain, binDmc);
		case 'bisection'
			threshold = threshold_bisection(symbolRatio, receivePower, equivalentDistribution, thresholdDomain, binDmc);
		case 'ml'
			threshold = threshold_ml(symbolRatio, receivePower);
		end

		% * Construct DMTC and recover i/o mapping
		dmac = dmc_integration(symbolRatio, receivePower, threshold);
		[~, sortIndex] = sort(receivePower);
		dmac(:, sortIndex) = dmac;

		% * Test convergence
		[wsr, rate] = rate_weighted(weight, snr, equivalentDistribution, dmac);
		isConverged = (wsr - wsr_) / wsr <= tolerance || isnan(wsr);
	end

	% * Update initializer
	Initializer.distribution = distribution;
	Initializer.beamforming = beamforming;
end
