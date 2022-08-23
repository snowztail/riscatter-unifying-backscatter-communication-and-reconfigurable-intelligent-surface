function [beamforming] = beamforming_dmrt(transmitPower, directChannel)
	% Function:
	%	- obtain MRT beamformer to direct channel
    %
    % Input:
	%	- transmitPower: average transmit power
	%	- directChannel [nTxs x 1]: direct AP-user channel
    %
    % Output:
	%	- beamforming [nTxs x 1]: transmit beamforming vector
    %
	% Comment:
	%	- optimal when no tags are present
	%
    % Author & Date: Yang (i@snowztail.com), 22 Jul 25

	beamforming = sqrt(transmitPower) * directChannel / norm(directChannel);
end
