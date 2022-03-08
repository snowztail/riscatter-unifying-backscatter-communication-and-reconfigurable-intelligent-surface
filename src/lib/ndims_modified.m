function [nDimensions] = ndims_modified(array)
	% Function:
    %	- get the actual dimensions of an array
	%	- obtain the optimal tag input distribution for a given discrete memoryless MAC
    %
    % Input:
	%	- array: the input array
    %
    % Output:
	%	- nDimensions: the actual dimensions of the array
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
