function [distribution, equivalentDistribution, convergence] = distribution_kkt(nTags, weight, snr, dmac, tolerance)
	% Function:
	%	- obtain KKT tag input probability distribution for maximizing weighted sum of primary and total backscatter rate
    %
    % Input:
	%	- nTags: number of tags
	%	- weight: relative priority of primary link
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
    %	- tolerance: minimum rate gain ratio per iteration
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- convergence [1 x nIterations]: weighted sum-rate at each iteration
    %
    % Comment:
    %	- iteratively update input distribution of each tag by coordinate descent
	%	- converge to KKT solution (stationary point)
	%	- DMAC is in joint form P(y | x_1, ..., x_K) instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Set default tolerance
	arguments
		nTags;
		weight;
		snr;
		dmac;
		tolerance = 1e-9;
	end

	% * Get data
	nInputs = size(dmac, 1);
	nStates = nthroot(nInputs, nTags);

	% ! Initialize distribution by previous solution
	persistent Initializer

	% * No previous solution, use uniform initializer
	if isempty(Initializer)
		Initializer.distribution = distribution_uniform(nTags, nStates);
	end

	% * Apply initializer
	distribution = Initializer.distribution;
	distributionTuple = tuple_tag(distribution);
	equivalentDistribution = prod(distributionTuple, 2);
	informationFunction = weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
	marginalInformation = marginal_information(distributionTuple, informationFunction);
	wsr = equivalentDistribution' * informationFunction;
	convergence = wsr;

	% * Iteratively update input distribution for all tags
	isConverged = false;
	while ~isConverged
		wsr_ = wsr;
		% * Update input distribution, information functions associated with each codeword, marginal information of each codeword, and mutual information for each tag
		for iTag = 1 : nTags
			distribution(:, iTag) = distribution(:, iTag) .* exp(marginalInformation(:, iTag)) ./ (distribution(:, iTag)' * exp(marginalInformation(:, iTag)));
			distributionTuple = tuple_tag(distribution);
			equivalentDistribution = prod(distributionTuple, 2);
			informationFunction = weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
			marginalInformation = marginal_information(distributionTuple, informationFunction);
			wsr = equivalentDistribution' * informationFunction;
		end
		convergence = [convergence, wsr];
		isConverged = (wsr - wsr_) / wsr <= tolerance || isnan(wsr);
	end

	% * Update initializer
	Initializer.distribution = distribution;
end


function [marginalInformation] = marginal_information(distributionTuple, informationFunction)
	% * Get data
	[nInputs, nTags] = size(distributionTuple);
	nStates = nthroot(nInputs, nTags);

	% * Construct tag state tuple
	stateTuple = tuple_tag(repmat(transpose(1 : nStates), [1, nTags]));

	% * Compute marginal information of each state of each tag
	marginalInformation = zeros(nStates, nTags, nInputs);
	for iState = 1 : nStates
		for iTag = 1 : nTags
			iInput = find(stateTuple(:, iTag) == iState);
			marginalInformation(iState, iTag, iInput) = prod(distributionTuple(iInput, [1 : iTag - 1, iTag + 1 : nTags]), 2) .* informationFunction(iInput);
		end
	end
	marginalInformation = sum(marginalInformation, 3);
end
