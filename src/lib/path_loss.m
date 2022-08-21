function [pathLoss] = path_loss(distance, exponent, rpl)
	% Function:
	%	- compute the path loss between transmitter and receiver
    %
    % Input:
    %	- distance: distance between transmitter and receiver
	%	- exponent: path loss exponent
	%	- rpl: reference path loss at 1m
    %
    % Output:
	%	- pathLoss: large-scale signal attenuation
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Set reference path loss at 1m
	arguments
		distance;
		exponent;
		rpl = db2pow(-30);
	end

	% * Compute path loss
	pathLoss = rpl * (1 / distance) ^ exponent;
end
