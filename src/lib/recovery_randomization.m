function [distribution, equivalentDistribution] = recovery_randomization(weight, snr, jointDistribution, dmac, nKernels, nInstances, tolerance)
	% Function:
	%	- extract individual tag input distribution from joint input distribution by randomization
    %
    % Input:
	%	- weight: relative priority of primary link
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- jointDistribution [nStates x ... (nTags-dimensional) ... x nStates]: joint tag input distribution with full transmit cooperation
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
	%	- nKernels: number of kernels used in appriximate array (i.e. rank of approximate array)
	%	- nInstances: number of random instances to generate for each tag
	%	- tolerance: minimum reduce of tolerable relative entropy (between original and approximate arrays) per iteration
    %
    % Output:
	%	- distribution [nStates x nTags]: tag input (i.e., state) probability distribution
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
    %
    % Comment:
	%	- approximate any joint input distribution array by rank-1 distribution in the mean sense
	%	- first approximate joint input distirbution by convex combination of rank-one kernels under minimum divergence criterion
    %	- then, for each kernel vector, generate random probability vectors whose mean equal to kernel vector
    %	- then, generate random probability vector groups whose mean equal to kernel vectors
	%	- select the probability vector group that maximizes weighted sum-rate
	%	- equivalent to a random distribution search with guiduance on correlation matrix
	%	- number of kernels and random samples are designable with performance-complexity tradeoff
	%	- candidates follows uniform distribution within a sphere bounded by probability simplex
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 03


	% * Declare default tolerance
	arguments
		weight;
		snr;
		jointDistribution;
		dmac;
		nKernels = 10;
		nInstances = 5e2;
		tolerance = 1e-9;
	end

	% * Get data
	nTags = dimension_number(jointDistribution);
	nStates = size(jointDistribution, 1);

	% * Initialize kernel vectors and relative entropy
	kernelVector = normalize(ones(nStates, nTags, nKernels), 'norm', 1);
	relativeEntropy = 1;

	% * Block coordinate descent (kernel coefficients -> kernel vectors of tag 1 -> kernel vectors of tag 2 -> ...)
	isConverged = false;
	iter = 0;
	while ~isConverged
		% * Update iteration index
		iter = iter + 1;
		relativeEntropy_ = relativeEntropy;

		% * Update kernel coefficients
		cvx_begin
			variable kernelCoefficient(nKernels, 1);
			expression approximateDistribution(size(jointDistribution));
			for iKernel = 1 : nKernels
				approximateDistribution = approximateDistribution + kernelCoefficient(iKernel) * product_outer(kernelVector(:, :, iKernel));
			end
			relativeEntropy = sum(vec(rel_entr(jointDistribution, approximateDistribution)));
			minimize relativeEntropy
			subject to
				kernelCoefficient == simplex(nKernels);
		cvx_end

		% * Update kernel vectors
		for iTag = 1 : nTags
			cvx_begin
				variable newVector(nStates, nKernels);
				expression approximateDistribution(size(jointDistribution));
				kernelVector = cvx(kernelVector);
				kernelVector(:, iTag, :) = newVector;
				for iKernel = 1 : nKernels
					approximateDistribution = approximateDistribution + kernelCoefficient(iKernel) * product_outer(kernelVector(:, :, iKernel));
				end
				relativeEntropy = sum(vec(rel_entr(jointDistribution, approximateDistribution)));
				minimize relativeEntropy
				subject to
					for iKernel = 1 : nKernels
						newVector(:, iKernel) == simplex(nStates);
					end
			cvx_end
		end

		% * Check convergence
		isConverged = (relativeEntropy - relativeEntropy_) <= tolerance || iter >= 1e2;
	end

	% * Generate random probability vectors with mean equal to kernel vectors
	wsr = 0;
	for iKernel = 1 : nKernels
		for iInstance = 1 : kernelCoefficient(iKernel) * nInstances
			distributionInstance = zeros(nStates, nTags);
			for iTag = 1 : nTags
				meanVector = kernelVector(:, iTag, iKernel);
				radius = min([abs(sum(meanVector) - meanVector - 1) / sqrt(nStates - 1); norm(meanVector) ^ 2]);

				% * Draw random vector uniformly within sphere
				isValid = false;
				while ~isValid
					randomVector = rand(nStates, 1) * 2 * radius - radius + meanVector;
					isValid = norm(randomVector - meanVector) <= radius;
				end
				distributionInstance(:, iTag) = randomVector + (1 - ones(1, nStates) * randomVector) / nStates * ones(nStates, 1);
			end
			equivalentDistributionInstance = prod(tuple_tag(distributionInstance), 2);
			wsrInstance = equivalentDistributionInstance' * (weight * information_primary(snr) + (1 - weight) * information_backscatter(equivalentDistributionInstance, dmac));
			if wsrInstance > wsr
				distribution = distributionInstance;
				equivalentDistribution = equivalentDistributionInstance;
				wsr = wsrInstance;
			end
		end
	end
end


function [outerProduct] = product_outer(matrix)
	nVectors = size(matrix, 2);
    outerProduct = matrix(:, 1);
	for iVector = 2 : nVectors
		vector = permute(matrix(:, iVector), circshift(1 : (dimension_number(outerProduct) + 1), dimension_number(outerProduct)));
		outerProduct = repmat(outerProduct, size(vector)) .* repmat(vector, size(outerProduct));
	end
end
