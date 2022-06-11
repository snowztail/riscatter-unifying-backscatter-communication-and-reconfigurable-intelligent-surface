function [distribution, equivalentDistribution] = distribution_exhaustion(nTags, weight, snr, dmac, resolution)
	% Function:
	%	- exhaustively search tag input probability distribution for maximizing weighted sum of primary and total backscatter rate
    %
    % Input:
	%	- nTags: number of tags
	%	- weight: relative priority of primary link
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
	%	- resolution: gap between neighbor search levels
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Set default resolution
	arguments
		nTags;
		weight;
		snr;
		dmac;
		resolution = 5e-3;
	end

	% * Get data
	nInputs = size(dmac, 1);
	nStates = nthroot(nInputs, nTags);

	% * Obtain probability simplex and tag-state tuples to evaluate
	[probabilitySimplex, nVertices] = simplex_probability(nStates, resolution);
	candidateTuple = tuple_tag(repmat(transpose(1 : nVertices), [1, nTags]));
	nTuples = size(candidateTuple, 1);

	% * Exhaustive search
	wsr = 0;
	for iTuple = 1 : nTuples
		distributionInstance = probabilitySimplex(:, candidateTuple(iTuple, :));
		equivalentDistributionInstance = prod(tuple_tag(distributionInstance), 2);
		wsrInstance = equivalentDistributionInstance' * (weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistributionInstance, dmac));
		if wsrInstance > wsr
			distribution = distributionInstance;
			equivalentDistribution = equivalentDistributionInstance;
			wsr = wsrInstance;
		end
	end
end


function [probabilitySimplex, nVertices] = simplex_probability(nStates, resolution)
	nBins = ceil(1 / resolution);
	stateTuple = transpose(nchoosek(1 : nStates + nBins - 1, nBins) - (0 : nBins - 1));
	nVertices = size(stateTuple, 2);
	probabilitySimplex = zeros(nStates, nVertices);
	for iState = 1 : nStates
		probabilitySimplex(iState, :) = sum(stateTuple == iState) / nBins;
	end
end
