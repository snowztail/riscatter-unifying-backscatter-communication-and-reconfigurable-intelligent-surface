function [nDimensions] = dimension_number(array)
	% Function:
    %	- return number of dimensions of an array
    %
    % Input:
	%	- array: input array
    %
    % Output:
	%	- nDimensions: number of dimensions of an array
    %
    % Comment:
	%	- return 1 for row/column vectors and 0 for []
    %	- https://stackoverflow.com/a/53694167/8891100
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 02

	if numel(array) == 1
		nDimensions = 1;
	else
		nDimensions = sum(size(array) > 1);
	end
end
