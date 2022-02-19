function [constellation] = psk(nStates, offset)
	% Function:
    %   - simulate the constellation diagram of M-ary phase shift keying
    %
    % Input:
    %   - nStates: number of constellation points
    %
    % Output:
	%	- constellation [nStates * 1]: the normalized coordinate of constellation diagram
	%	- offset: the phase offset
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

	% * Set default offset
	arguments
		nStates;
		offset = 0;
	end

	% * Simulate the constellation diagram
	constellation = exp(1i * (2 * pi * (1 : nStates)) / nStates + offset);
end
