function [jointDistribution, equivalentDistribution] = distribution_cooperation(nTags, weight, snr, dmac, tolerance)
	% Function:
	%	- obtain the optimal joint tag input distribution with full transmit cooperation for maximizing weighted sum of primary and total backscatter rate
    %
    % Input:
	%	- nTags: number of tags
	%	- weight: relative priority of primary link
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
    %	- tolerance: minimum rate gain ratio per iteration
    %
    % Output:
	%	- jointDistribution [nStates x ... (nTags-dimensional) ... x nStates]: joint tag input distribution with full transmit cooperation
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Comment:
	%	- with full transmit cooperation, the tags are equivalent to a single information source with input distribution in nInputs-dimensional probability simplex
    %   - in general, tags encode independently and it is hard to approximate the joint distribution array without cooperation
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 17

	% * Set default tolerance
	arguments
		nTags;
		weight;
		snr;
		dmac;
		tolerance = 1e-6;
	end

	% * Get data
	nInputs = size(dmac, 1);
	nStates = nthroot(nInputs, nTags);

	% ! Initialize distribution by previous solution
	persistent Initializer

	% * No previous solution, use uniform initializer
	if isempty(Initializer)
		Initializer.equivalentDistribution = normalize(ones(nStates ^ nTags, 1), 'norm', 1);
	end

	% * Apply initializer
	equivalentDistribution = Initializer.equivalentDistribution;
	informationFunction = weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
	wsr = equivalentDistribution' * informationFunction;

	% * Iteratively update equivalent distribution, information function associated with each (joint) codeword, and weighted sum achievable rate
	isConverged = false;
	while ~isConverged
		wsr_ = wsr;
		equivalentDistribution = equivalentDistribution .* exp(informationFunction) / (equivalentDistribution' * exp(informationFunction));
		jointDistribution = permute(reshape(equivalentDistribution, nStates * ones(1, nTags)), nTags : -1 : 1);
		informationFunction = weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
		wsr = equivalentDistribution' * informationFunction;
		isConverged = (wsr - wsr_) / wsr <= tolerance || isnan(wsr);
	end

	% * Update initializer
	Initializer.equivalentDistribution = equivalentDistribution;
end
