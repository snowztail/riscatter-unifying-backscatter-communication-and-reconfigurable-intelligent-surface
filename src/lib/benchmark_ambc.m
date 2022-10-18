function [rate, distribution, threshold] = benchmark_ambc(nTags, symbolRatio, transmitPower, noisePower, equivalentChannel, scatterRatio, directChannel, cascadedChannel)
	% Function:
	%	- evaluate the achievable primary and total backscatter rate for ambient backscatter communications
    %
    % Input:
	%	- nTags: number of tags
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
	%	- scatterRatio: amplitude scatter ratio at tags
	%	- directChannel [nTxs x 1]: AP-user channel
	%	- cascadedChannel [nTxs x nTags]: cascaded forward-backward channel of all tags
    %
    % Output:
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
    % Comment:
	%	- assume single Tx, uniform input distribution and maximum likelihood energy detection (no need for BCD)
	%	- consider co-located non-cooperative receivers as benchmark
	%	- use SINR to approximate primary rate even if backscatter sources are not Gaussian and scattered signals are correlated with primary source
    %
    % Author & Date: Yang (i@snowztail.com), 22 Oct 14

	% * Get data
	nInputs = size(equivalentChannel, 2);
	nStates = nthroot(nInputs, nTags);

	% * No beamformer and equiprobable inputs
	distribution = distribution_uniform(nTags, nStates);
	beamforming = sqrt(transmitPower);
	equivalentDistribution = prod(tuple_tag(distribution), 2);

	% * Compute average receive power and SINR
	receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
	sinr = abs(directChannel' * beamforming) .^ 2 / sum(abs(scatterRatio * cascadedChannel' * beamforming) .^ 2 + noisePower) * ones(nInputs, 1);

	% * Obtain ML decision threshold
	threshold = threshold_ml(symbolRatio, receivePower);

	% * Construct DMTC and recover i/o mapping
	dmac = dmc_integration(symbolRatio, receivePower, threshold);
	[~, sortIndex] = sort(receivePower);
	dmac(:, sortIndex) = dmac;

	% * Compute achievable rates
	[~, rate] = rate_weighted(sinr, equivalentDistribution, dmac);
end
