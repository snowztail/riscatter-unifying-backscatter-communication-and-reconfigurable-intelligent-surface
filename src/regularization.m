function [regularizedArray] = regularization(array)
	% Function:
	%	- discard the negative elements and normalize the array to unit sum
    %
    % Input:
	%	- array: the input array
    %
    % Output:
	%	- regularizedArray: regularized array with unit sum and no negative elements
    %
    % Comment:
    %	- deal with the precision issue
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 10

	array(array < 0) = 0;
	regularizedArray = array ./ sum(vec(array));
end
