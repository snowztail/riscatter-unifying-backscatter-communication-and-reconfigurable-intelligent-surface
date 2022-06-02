function [directDistance, forwardDistance, backwardDistance] = layout(nTags)
	% Function:
	%	- compute distance of direct, forward, and backward paths
    %
    % Input:
	%	- nTags: number of tags
    %
    % Output:
	%	- directDistance: distance of direct path
	%	- forwardDistance: distance of forward path
	%	- backwardDistance: distance of backward path
    %
	% Comment:
	%	- layout is customizable
	%
    % Author & Date: Yang (i@snowztail.com), 22 Jun 02

	% * Parameters
	coverage = 10;
	radius = 1;

	% * Coordinates
	ap = [coverage; 0];
	user = [0; 0];
	[tag(1, :), tag(2, :)] = pol2cart(2 * pi * rand(1, nTags), sqrt(radius * rand(1, nTags)));
	% tag(1, :) = ones(1, nTags);
	% tag(2, :) = (1 : nTags) - mean(1 : nTags);

	% * Distances
	directDistance = norm(ap - user);
	forwardDistance = vecnorm(ap - tag);
	backwardDistance = vecnorm(tag - user);
end
