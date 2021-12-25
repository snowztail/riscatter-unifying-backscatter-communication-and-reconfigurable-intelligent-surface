function [combination] = nested_combvec(vector, times)
	% Function:
    %   - create all combinations of copies of a vector
    %
    % Input:
    %   - vector: the vector to be repeated and combined (as the argument in nested for loops)
	%	- times: the number of copies to combine
    %
    % Output:
	%	- combination: the combination matrix whose columns consist of all combinations
    %
    % Author & Date: Yang (i@snowztail.com), 21 Dec 25

	if times > 1
		combination = sortrows(combvec(vector, nested_combvec(vector, times - 1)));
	elseif times == 1
		combination = vector;
	end
end
