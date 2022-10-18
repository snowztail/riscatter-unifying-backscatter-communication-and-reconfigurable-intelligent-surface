function [rate, distribution] = benchmark_sr(nTags, transmitPower, noisePower, equivalentChannel)
	% Function:
	%	- evaluate the achievable primary and total backscatter rate for commensal symbiotic ratio
    %
    % Input:
	%	- nTags: number of tags
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
    %
    % Output:
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
    %
    % Comment:
	%	- assume single Tx, uniform input distribution, and ideal primary detection, SIC, backscatter detection from intermediate signal (no need for BCD)
	%	- when symbol radio is sufficiently large:
	%		- the primary ergodic rate under semi-coherent detection can be approximated by that of coherent detection
	%		- the total backscatter achievable rate approaches noise-free case
    %
    % Author & Date: Yang (i@snowztail.com), 22 Oct 17

	% * Get data
	nInputs = size(equivalentChannel, 2);
	nStates = nthroot(nInputs, nTags);

	% * No beamformer and equiprobable inputs
	distribution = distribution_uniform(nTags, nStates);
	beamforming = sqrt(transmitPower);
	equivalentDistribution = prod(tuple_tag(distribution), 2);

	% * Compute signal power and SNR
	signalPower = abs(equivalentChannel' * beamforming) .^ 2;
	snr = signalPower / noisePower;

	% * Compute asymptotical achievable rates
	primaryRate = equivalentDistribution' * log(1 + snr);
	backscatterRate = nTags * log(nStates);
	rate = [primaryRate; backscatterRate];
end
