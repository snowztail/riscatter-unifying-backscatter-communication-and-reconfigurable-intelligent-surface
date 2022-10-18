function [rate, distribution] = benchmark_ris(nTags, transmitPower, noisePower, equivalentChannel)
	% Function:
	%	- evaluate the achievable primary and total backscatter rate for reconfigurable intelligent surface
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
	%	- assume single Tx and degenrate input distribution
	%	- primary equivalent channel is constant within each block
    %
    % Author & Date: Yang (i@snowztail.com), 22 Oct 17

	% * Get data
	nInputs = size(equivalentChannel, 2);
	nStates = nthroot(nInputs, nTags);

	% * No beamformer
	beamforming = sqrt(transmitPower);

	% * Compute signal power and SNR
	signalPower = abs(equivalentChannel' * beamforming) .^ 2;
	snr = signalPower / noisePower;

	% * Determine RIS state
	equivalentDistribution = (signalPower == max(signalPower));
	stateTuple = tuple_tag(repmat(transpose(1 : nStates), [1, nTags]));
	distribution = zeros(nStates, nTags);
	for iTag = 1 : nTags
		distribution(stateTuple(equivalentDistribution, iTag), iTag) = 1;
	end

	% * Compute achievable rates
	primaryRate = equivalentDistribution' * log(1 + snr);
	backscatterRate = 0;
	rate = [primaryRate; backscatterRate];
end
