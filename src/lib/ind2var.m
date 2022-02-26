function [variable] = ind2var(prefix, suffix)
	% Function:
	%	- combine prefix and suffix to create a string as variable name
    %
    % Input:
	%	- prefix: character or string as prefix of variable name
	%	- suffix: index as suffix of variable
    %
    % Output:
	%	- variable: variable name with prefix and suffix
    %
    % Comment:
    %   - tailored hack for CVX to create different variable names
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 26

	variable = strcat(prefix, sprintf('%d', suffix));
end
