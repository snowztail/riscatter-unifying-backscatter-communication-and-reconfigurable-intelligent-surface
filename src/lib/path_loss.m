function [pathLoss] = path_loss(frequency, distance, exponent, txGain, rxGain)
	% Function:
	%	- compute the path loss between transmitter and receiver
    %
    % Input:
	%	- frequency: carrier frequency
    %	- distance: distance between transmitter and receiver
	%	- exponent: path loss exponent
	%	- txGain: transmit antenna gain
	%	- rxGain: receive antenna gain
    %
    % Output:
	%	- pathLoss: large-scale signal attenuation
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Set default antenna gains
	arguments
		frequency;
		distance;
		exponent;
		txGain = 1;
		rxGain = 1;
	end

	% * Compute path loss
	wavelength = physconst('LightSpeed') / frequency;
	pathLoss = txGain * rxGain * wavelength ^ 2 / (16 * pi ^ 2 * distance ^ exponent);
end
