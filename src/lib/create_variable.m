function [cvx_problem] = create_variable(name, size, cvx_problem)
	% Function:
	%	- create CVX variable with specified name and size
    %
    % Input:
	%	- name: variable name
	%	- size: variable size
    %
    % Comment:
    %   - tailored hack for CVX to create variable in loops
	%
    % Author & Date: Yang (i@snowztail.com), 22 Feb 26

	eval(['variable', ' ', name, '(', size, ')']);
end
