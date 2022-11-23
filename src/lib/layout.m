function [forwardDistance, backwardDistance] = layout(directDistance, nTags, coverage)
	% Function:
	%	- compute distance of direct, forward, and backward paths
    %
    % Input:
	%	- directDistance: AP-user distance
	%	- nTags: number of tags
	%	- coverage [1 x nTags]: maximum tag(s)-user distance
    %
    % Output:
	%	- forwardDistance [1 x nTags]: AP-tag(s) distance
	%	- backwardDistance [1 x nTags]: tag(s)-user distance
    %
	% Comment:
	%	- assume tags are uniformly dropped within a circle centered at user
	%	- tags farther from AP/user contribute less to equivalent channel
	%
    % Author & Date: Yang (i@snowztail.com), 22 Jun 02

	% * Coordinates
	ap = [directDistance; 0];
	user = [0; 0];
	[tag(1, :), tag(2, :)] = pol2cart(2 * pi * rand(1, nTags), sqrt(coverage ^ 2 * rand(1, nTags)));
	% tag(1, :) = ones(1, nTags);
	% tag(2, :) = (1 : nTags) - mean(1 : nTags);

	% * Distances
	forwardDistance = vecnorm(ap - tag);
	backwardDistance = vecnorm(tag - user);
end
