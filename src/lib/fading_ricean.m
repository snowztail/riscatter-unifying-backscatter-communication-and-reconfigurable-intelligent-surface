function [fading] = fading_ricean(nTxs, nRxs, factor, aoa)
	% Function:
	%	- simulate small-scale fading in i.i.d. Ricean distribution
    %
    % Input:
	%	- nTxs: number of transmit antennas
	%	- nRxs: number of receive antennas
	%	- factor: Ricean factor
	%	- aoa: angle of arrival
    %
    % Output:
	%	- fading [nTxs x nRxs]: small-scale signal attenuation
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Set default angle of arrival
	arguments
		nTxs;
		nRxs;
		factor;
		aoa = rand * 2 * pi;
	end

	% * Simulate Ricean fading
	fading = sqrt(0.5 * factor / (factor + 1)) * exp(1i * aoa) + sqrt(0.5 / (factor + 1)) * (randn(nTxs, nRxs) + 1i * randn(nTxs, nRxs));
end
