function [instance] = rand_normalized(arraySize, dimension)
	% Function:
    %	- generate a random array with all elements in uniform distribution and normalized sum in specific dimension
    %
    % Input:
    %	- arraySize: size of the desired array
	%	- dimension: the dimension where all elements sum up to 1
    %
    % Output:
	%	- instance: an instance of the random array with all elements in uniform distribution and normalized sum in specific dimension
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 09

	instance = rand(arraySize);
	instance = instance ./ sum(instance, dimension);
end
