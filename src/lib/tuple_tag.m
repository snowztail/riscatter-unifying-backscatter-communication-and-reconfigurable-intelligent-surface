function [tuple] = tuple_tag(domain)
	% Function:
	%	- enumerate tag state index/input probability tuples
    %
    % Input:
	%	- domain [nStates x nTags]: k-th column denotes state index set or input probability set of tag k
    %
    % Output:
	%	- tuple [nInputs x nTags]: i-th row is a K-tuple whose entries come from state index/input distribution of all tags
    %
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	% * Get data
	nTags = size(domain, 2);

	% * Enumerate tuples by recursion
	if nTags > 1
		tuple = product_cartesian(domain(:, 1), tuple_tag(domain(:, 2 : end)));
	elseif nTags == 1
		tuple = domain;
	end
end
