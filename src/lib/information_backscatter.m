function [backscatterInformation] = information_backscatter(equivalentDistribution, dmac)
	% Function:
	%	- compute total backscatter information function for each tag state tuple
    %
    % Input:
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
    %
    % Output:
	%	- backscatterInformation [nInputs x 1]: total backscatter information function for each tag state tuple
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 19

	% * Get data
	[nInputs, nOutputs] = size(dmac);

	% * Compute total backscatter information function
	backscatterInformation = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformation(iInput, iOutput) = dmac(iInput, iOutput) * log(dmac(iInput, iOutput) / (equivalentDistribution' * dmac(:, iOutput)));
		end
	end
	backscatterInformation = sum(backscatterInformation, 2);
end
