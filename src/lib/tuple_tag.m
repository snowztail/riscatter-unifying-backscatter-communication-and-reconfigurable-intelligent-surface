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

function [tuple] = product_cartesian(varargin)
	% Function:
	%	- compute (N-fold) Cartesian product of N input sets
    %
    % Input:
	%	- varargin: each input set is a column vector
    %
    % Output:
	%	- tuple: each row is a N-tuple whose n-th entry come from set C_n
    %
	% Comment:
	%	- Cartesian product of sets C_1, ..., C_N is defined as the set of all ordered N-tuples (x_1, ..., x_N) such that x_n belongs to C_n
	%	- modified from built-in function `combvec`
	%
    % Author & Date: Yang (i@snowztail.com), 22 May 10

	if isempty(varargin)
		tuple = [];
	else
		tuple = varargin{1};
		for i = 2 : length(varargin)
			tail = varargin{i};
			% tuple = [repelem(tuple, size(tail, 1), size(tail, 2)), copy_interleaved(tail, size(tuple, 1))];
			tuple = [repelem(tuple, size(tail, 1), 1), repmat(tail, size(tuple, 1), 1)];
		end
	end
end
