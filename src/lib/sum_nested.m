function [marginalArray] = sum_nested(array, dimension)
	% Function:
    %	- obtain the sum of array over specific dimensions
    %
    % Input:
    %	- array: the input array
	%	- dimension: the dimensions to be summed
    %
    % Output:
	%	- marginalArray: the sum of array over specific dimensions
	%
	% Comment:
	%	- CVX does not support sum over multiple dimensions yet
	%	- Require dimension indexes in ascending order
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 01

	marginalArray = array;
	if length(dimension) >= 1
		marginalArray = sum_nested(sum(marginalArray, dimension(end)), dimension(1 : end - 1));
	end
end
