function [distribution, equivalentDistribution] = distribution_kkt(nTags, symbolRatio, weight, snr, dmac, tolerance)
	% Function:
	%	- obtain KKT tag input probability distribution for maximizing weighted sum of primary and total backscatter rate
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
    %	- iteratively update input distribution of each tag by coordinate descent
	%	- converge to KKT solution (stationary point)
	%	- DMAC is in joint form P(y | x_1, ..., x_K) instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Jan 09

	% * Set default tolerance
	arguments
		nTags;
		symbolRatio;
		weight;
		snr;
		dmac;
		tolerance = 1e-10;
	end

	% * Get data
	nInputs = size(dmac, 1);
	nStates = nthroot(nInputs, nTags);

	% * Initialize input distribution as uniform
	distribution = normalize(ones(nStates, nTags), 'norm', 1);
	distributionTuple = tuple_tag(distribution);
	equivalentDistribution = prod(distributionTuple, 2);
	informationFunction = weight * information_primary(symbolRatio, snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
	marginalInformation = marginal_information(distributionTuple, informationFunction);
	wsr = equivalentDistribution' * informationFunction;

	% * Iteratively update input distribution for all tags
	isConverged = false;
	while ~isConverged
		wsr_ = wsr;
		% * Update input distribution, information functions associated with each codeword, marginal information of each codeword, and mutual information for each tag
		for iTag = 1 : nTags
			distribution(:, iTag) = distribution(:, iTag) .* exp(marginalInformation(:, iTag)) ./ (distribution(:, iTag)' * exp(marginalInformation(:, iTag)));
			distributionTuple = tuple_tag(distribution);
			equivalentDistribution = prod(distributionTuple, 2);
			informationFunction = weight * information_primary(symbolRatio, snr) + (1 - weight) * information_backscatter(equivalentDistribution, dmac);
			marginalInformation = marginal_information(distributionTuple, informationFunction);
			wsr = equivalentDistribution' * informationFunction;
		end
		isConverged = abs(wsr - wsr_) <= tolerance;
	end
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
			marginalInformation(iState, iTag, iInput) = prod(distributionTuple(iInput, setdiff(1 : nTags, iTag)), 2) .* informationFunction(iInput);
		end
	end
	marginalInformation = sum(marginalInformation, 3);
end
