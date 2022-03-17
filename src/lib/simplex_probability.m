function [probabilitySimplex] = simplex_probability(nStates, resolution)
	% Function:
    %	- enumerate all probability distribution candidates under specified resolution
    %
    % Input:
    %	- nStates: number of states in tag constellation diagram
	%	- resolution: minimum gap between quantization levels
    %
    % Output:
	%	- probabilitySimplex [nCandidates * nStates]: each row represents a candidate probability distribution
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Declare default resolution
	arguments
		nStates;
		resolution = 1e-2;
	end

	% * Obtain the number of discretization bins
	nBins = 1 / resolution;

	% * Label the state to which each region contributes
	label = nchoosek(1 : nStates + nBins - 1, nBins) - (0 : nBins - 1);

	% * Evaluate the probality of occurrence
	probabilitySimplex = zeros(size(label, 1), nStates);
	for iState = 1 : nStates
		probabilitySimplex(:, iState) = flipud(sum(label == iState, 2) / nBins);
	end
end
