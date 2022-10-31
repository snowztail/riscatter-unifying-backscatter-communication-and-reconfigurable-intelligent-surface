function [rate] = benchmark_legacy(transmitPower, noisePower, directChannel)
	% Function:
	%	- evaluate the achievable primary and total backscatter rate for legacy transmission
    %
    % Input:
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- directChannel [nTxs x 1]: direct AP-user channel
    %
    % Output:
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
    %
    %
    % Author & Date: Yang (i@snowztail.com), 22 Oct 31

	% * No beamformer
	beamforming = sqrt(transmitPower);

	% * Compute signal power and SNR
	signalPower = abs(directChannel' * beamforming) .^ 2;
	snr = signalPower / noisePower;

	% * Compute achievable rates
	primaryRate = log(1 + snr);
	backscatterRate = 0;
	rate = [primaryRate; backscatterRate];
end
