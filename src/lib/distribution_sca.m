function [distribution, equivalentDistribution] = distribution_sca(nTags, symbolRatio, weight, snr, dmac, tolerance)
	% Function:
	%	- optimize tag input probability distribution by successive convex approximation for maximizing weighted sum of primary and total backscatter rate
    %
    % Input:
	%	- nTags: number of tags
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- weight: relative priority of primary link
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Comment:
	%	- sequentially approximate the equivalent (i.e. product of individual) distribution by its first-order Taylor expansion at previous point
	%	- expected to converge to KKT points but convergence is unstable
    %
    % Author & Date: Yang (i@snowztail.com), 23 Mar 15

	% * Set default tolerance
	arguments
		nTags;
		symbolRatio;
		weight;
		snr;
		dmac;
		tolerance = 1e-9;
	end

	% * Get data
	nInputs = size(dmac, 1);
	nStates = nthroot(nInputs, nTags);

	% * Construct tag state tuple
	stateTuple = tuple_tag(repmat(transpose(1 : nStates), [1, nTags]));

	% * Initialize input distribution as uniform
	distribution = normalize(ones(nStates, nTags), 'norm', 1);
	equivalentDistribution = prod(tuple_tag(distribution), 2);
	wsr = rate_weighted(symbolRatio, weight, snr, equivalentDistribution, dmac);

	% * Iteratively update input distribution by SCA
	isConverged = false;
	while ~isConverged
		distribution_ = distribution;
		wsr_ = wsr;
		cvx_begin
			variable distribution(nStates, nTags);
			expression equivalentDistributionSca(nInputs, 1);
			for iInput = 1 : nInputs
				for iTag = 1 : nTags
					equivalentDistributionSca(iInput) = equivalentDistributionSca(iInput) + distribution(stateTuple(iInput, iTag), iTag) * prod(distribution_(sub2ind(size(distribution_), stateTuple(iInput, setdiff(1 : nTags, iTag)), setdiff(1 : nTags, iTag))), 2);
				end
				equivalentDistributionSca(iInput) = equivalentDistributionSca(iInput) - (nTags - 1) * prod(distribution_(sub2ind(size(distribution_), stateTuple(iInput, :), 1 : nTags)), 2);
			end
			wsrSca = rate_weighted(symbolRatio, weight, snr, equivalentDistributionSca, dmac);
			maximize wsrSca
			subject to
				for iTag = 1 : nTags
					distribution(:, iTag) == simplex(nStates);
				end
		cvx_end

		% * Compute actual weighted sum rate
		equivalentDistribution = prod(tuple_tag(distribution), 2);
		wsr = rate_weighted(symbolRatio, weight, snr, equivalentDistribution, dmac);

		% * Test convergence
		isConverged = abs(wsr - wsr_) <= tolerance;
	end
end
