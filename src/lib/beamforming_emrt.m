function [beamforming] = beamforming_emrt(transmitPower, equivalentChannel, equivalentDistribution)
	% Function:
	%	- obtain MRT beamformer to equivalent primary channel
    %
    % Input:
	%	- transmitPower: average transmit power
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Output:
	%	- beamforming [nTxs x 1]: transmit beamforming vector
    %
	% Comment:
	%	- a balanced low-complexity scheme
	%
    % Author & Date: Yang (i@snowztail.com), 22 Jul 25

	ric = equivalentChannel * equivalentDistribution;
	beamforming = sqrt(transmitPower) * ric / norm(ric);
end
