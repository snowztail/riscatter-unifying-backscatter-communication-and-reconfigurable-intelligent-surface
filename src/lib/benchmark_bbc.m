function [rate, distribution, threshold] = benchmark_bbc(nTags, symbolRatio, transmitPower, noisePower, equivalentChannel)
	% Function:
	%	- evaluate the achievable primary and total backscatter rate for bistatic backscatter communications
    %
    % Input:
	%	- nTags: number of tags
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
    %
    % Output:
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
    % Comment:
	%	- assume single Tx, uniform input distribution and maximum likelihood energy detection (no need for BCD)
	%	- deterministic CW component can be either cancelled or modeled within translated constellation
    %
    % Author & Date: Yang (i@snowztail.com), 22 Oct 15

	% * Get data
	nInputs = size(equivalentChannel, 2);
	nStates = nthroot(nInputs, nTags);

	% * No beamformer and equiprobable inputs
	distribution = distribution_uniform(nTags, nStates);
	beamforming = sqrt(transmitPower);
	equivalentDistribution = prod(tuple_tag(distribution), 2);

	% * Compute signal power
	signalPower = abs(equivalentChannel' * beamforming) .^ 2;

	% * No primary signal, null primary SNR
	snr = zeros(nInputs, 1);

	% * Obtain ML decision threshold
	threshold = threshold_ml_bbc(symbolRatio, signalPower, noisePower);

	% * Construct DMTC and recover i/o mapping
	dmac = dmc_integration_bbc(symbolRatio, signalPower, noisePower, threshold);
	[~, sortIndex] = sort(signalPower);
	dmac(:, sortIndex) = dmac;

	% * Compute achievable rates
	[~, rate] = rate_weighted(snr, equivalentDistribution, dmac);
end

function [threshold] = threshold_ml_bbc(symbolRatio, signalPower, noisePower)
	nOutputs = size(signalPower, 1);

	sortPower = symbolRatio * sort(signalPower);
	threshold(1) = sortPower(1);
	threshold(nOutputs + 1) = icdf('Gamma', 1 - eps, symbolRatio, noisePower) + sortPower(end);
	for iOutput = 2 : nOutputs
		threshold(iOutput) = (sortPower(iOutput - 1) * exp((sortPower(iOutput - 1) - sortPower(iOutput)) / (noisePower * (symbolRatio - 1))) - sortPower(iOutput)) / (exp((sortPower(iOutput - 1) - sortPower(iOutput)) / (noisePower * (symbolRatio - 1))) - 1);
	end
end

function [dmc] = dmc_integration_bbc(symbolRatio, signalPower, noisePower, threshold)
	nInputs = size(signalPower, 1);
	nOutputs = size(threshold, 2) - 1;

	dmc = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		dmc(iInput, :) = diff(cdf('Gamma', threshold - symbolRatio * signalPower(iInput), symbolRatio, noisePower), [], 2);
	end
	dmc(dmc < eps) = eps;
end
