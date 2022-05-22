function [fading] = fading_rayleigh(nTxs, nRxs)
	% Function:
	%	- simulate small-scale fading in i.i.d. Rayleigh (CSCG) distribution with unit variance
    %
    % Input:
	%	- nTxs: number of transmit antennas
	%	- nRxs: number of receive antennas
    %
    % Output:
	%	- fading [nTxs x nRxs]: small-scale signal attenuation
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Simulate Rayleigh fading
	fading = sqrt(0.5) * (randn(nTxs, nRxs) + 1i * randn(nTxs, nRxs));
end
