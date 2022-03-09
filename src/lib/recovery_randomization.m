function [inputDistribution, equivalentDistribution, weightedSumRate] = recovery_randomization(jointDistribution, dmtc, weight, symbolRatio, snr, nKernels, nSamples, tolerance)
	% Function:
	%	- extract a good tag input distribution (corresponding to no transmit cooperation) by randomization
    %
    % Input:
	%	- jointDistribution [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
	%	- nKernels: number of kernels used in the approximation (equals to the rank)
	%	- nSamples: number of random samples to be generated in randomization
	%	- tolerance: the maximum tolerable relative entropy between the original and approximated arrays
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
	%	- approximate any joint input distribution array by reduced-rank distribution
    %	- generate random samples of candidate tag probability vectors whose outer product expectation equals to the joint input distribution array
	%	- equivalent to a random search with guiduance on the correlation matrix of the distributions
	%	- the number of kernels and random samples are designable (performance-complexity tradeoff)
	%	- the candidates follows uniform distribution within a sphere bounded within the probability simplex
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 3

	% * Declare default tolerance
	arguments
		jointDistribution;
		dmtc;
		weight;
		symbolRatio;
		snr;
		nKernels = 4;
		nSamples = 5e2;
		tolerance = 1e-8;
	end

	% * Get data
	nStates = size(jointDistribution, 1);
	nTags = ndims_modified(jointDistribution);

	% * Block coordinate descent
	isConverged = false;
	kernelBasis = rand_normalized([nStates, nTags, nKernels], 1);
	while ~isConverged
		% * Update kernel coefficients
		cvx_begin
			variable kernelCoefficient(nKernels, 1)
			referenceDistribution = zeros(size(jointDistribution));
			for iKernel = 1 : nKernels
				referenceDistribution = referenceDistribution + kernelCoefficient(iKernel) * outer_product(kernelBasis(:, :, iKernel));
			end
			relativeEntropy = sum(vec(rel_entr(jointDistribution, referenceDistribution)));
			minimize relativeEntropy
			subject to
				kernelCoefficient == simplex(nKernels);
		cvx_end

		% * Update kernel bases
		for iTag = 1 : nTags
			cvx_begin
				variable newBasis(nStates, nKernels)
				kernelBasis = cvx(kernelBasis);
				kernelBasis(:, iTag, :) = newBasis;
				referenceDistribution = zeros(size(jointDistribution));
				for iKernel = 1 : nKernels
					referenceDistribution = referenceDistribution + kernelCoefficient(iKernel) * outer_product(kernelBasis(:, :, iKernel));
				end
				relativeEntropy = sum(vec(rel_entr(jointDistribution, referenceDistribution)));
				minimize relativeEntropy
				subject to
					for iKernel = 1 : nKernels
						newBasis(:, iKernel) == simplex(nStates);
					end
			cvx_end
		end
		isConverged = abs(relativeEntropy) <= tolerance;
	end

	% * Generate random vectors with kernel bases as prescribed mean
	weightedSumRate = 0;
	for iKernel = 1 : nKernels
		for iSample = 1 : round(kernelCoefficient(iKernel) * nSamples)
			inputDistributionCandidate = zeros(nTags, nStates);
			for iTag = 1 : nTags
				% * Retrieve the radius of the uniform random vector whose boundary lies within the probability simplex
				meanVector = kernelBasis(:, iTag, iKernel);
				radiusBound = zeros(nStates + 1, 1);
				for iState = 1 : nStates
					radiusBound(iState) = abs(sum(meanVector(setdiff(1 : nStates, iState))) - 1) / sqrt(nStates - 1);
				end
				radiusBound(nStates + 1) = norm(meanVector) ^ 2;
				radius = min(radiusBound);
				% * Generate random vector uniformly distributed within the nStates-dimensional sphere
				isValid = false;
				while ~isValid
					randomVector = 2 * radius * rand(nStates, 1) - radius + meanVector;
					isValid = norm(randomVector - meanVector) <= radius;
				end
				inputDistributionCandidate(iTag, :) = transpose(randomVector + (1 - ones(1, nStates) * randomVector) / nStates * ones(nStates, 1));
			end
			equivalentDistributionCandidate = prod(combination_distribution(inputDistributionCandidate), 1);
			weightedSumRateCandidate = weighted_sum_rate(weight, symbolRatio, snr, equivalentDistributionCandidate, dmtc);
			if weightedSumRateCandidate > weightedSumRate
				inputDistribution = inputDistributionCandidate;
				equivalentDistribution = equivalentDistributionCandidate;
				weightedSumRate = weightedSumRateCandidate;
			end
		end
	end
end
