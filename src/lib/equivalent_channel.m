function [equivalentChannel] = equivalent_channel(directChannel, cascadedChannel, reflectRatio, constellation)
	% Function:
    %   - obtain the equivalent channel candidates for primary transmission
    %
    % Input:
    %   - directChannel [nTxs * nRxs]: the AP-UE channel
	%	- cascadedChannel [nTxs * nRxs]: the cascaded AP-TG and TG-UE channel without tag reflection
	%	- reflectRatio: power reflection efficiency of the tag at a given direction
	%	- constellation [nStates * 1]: the normalized coordinate of constellation diagram
    %
    % Output:
	%	- equivalentChannel [nStates * 1]: the equivalent channel candidates corresponding to tag modulation states
    %
    % Comment:
	%	- assume the AP, tag, and user are of single-antenna
	%	- different tag states lead to different equivalent channel
	%	- for AmBC, the maximum-likelihood estimator reduces to energy detector based on equivalent channel strength
    %   - once the tag is successfully decoded, the corresponding equivalent channel can be retrieved from the candidates
    %
    % Author & Date: Yang (i@snowztail.com), 21 Sep 14

	% * Obtain equivalent channel candidates corresponding to tag states
	equivalentChannel = directChannel + sqrt(reflectRatio) * cascadedChannel * constellation;
end
