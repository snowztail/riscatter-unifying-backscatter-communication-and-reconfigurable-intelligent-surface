function [powerSet, complementSet] = power_set(set)
	% Function:
    %   - obtain the power set (exclude the empty set) of the input set
    %
    % Input:
    %   - set [1 * nElements]: the input set as row vector
    %
    % Output:
	%	- powerSet: a cell array containing all subsets (exclude the empty set) of the input set
	%	- complementSet: a cell array containing the complement of all subsets (exclude the empty set) of the input set
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 25

	powerSet = [];
	for iSize = 0 : length(set)
		powerSet = [powerSet; num2cell(nchoosek(set, iSize), 2)];
	end
	complementSet = flipud(powerSet(1 : end - 1));
	powerSet = powerSet(2 : end);
end
