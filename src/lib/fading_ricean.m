function [fading] = fading_ricean(nTxs, nRxs, factor, phase)
	% Function:
	%	- simulate small-scale fading in i.i.d. Ricean distribution
    %
    % Input:
	%	- nTxs: number of transmit antennas
	%	- nRxs: number of receive antennas
	%	- factor: Ricean factor
	%	- phase: phase shift of direct path between transmitter and receiver
    %
    % Output:
	%	- fading [nTxs x nRxs]: small-scale signal attenuation
	%
	% Comment:
	%	- does not consider Doppler effect
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Set default angle of arrival
	arguments
		nTxs;
		nRxs;
		factor;
		phase = rand * 2 * pi;
	end

	% * Simulate Ricean fading
	fading = sqrt(factor / (factor + 1)) * exp(1i * phase) + sqrt(0.5 / (factor + 1)) * (randn(nTxs, nRxs) + 1i * randn(nTxs, nRxs));
end
