function [primaryInformation] = information_primary(snr)
	% Function:
	%	- compute primary information function for each tag state tuple
    %
    % Input:
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
    %
    % Output:
	%	- primaryInformation [nInputs x 1]: primary information function for each tag state tuple
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 19

	primaryInformation = log(1 + snr);
end
