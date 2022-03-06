function [randomizedInputDistribution, randomizedEquivalentDistribution, randomizedRate] = recovery_randomization(jointArray, dmtc, nKernels, nSamples, tolerance)
	% Function:
	%	- approximate any joint input distribution array by reduced-rank distribution
    %
    % Input:
	%	- jointArray [nStates * ... (nTags) ... * nStates]: the joint input distribution of all tags corresponding to the relaxed input optimization problem
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- nKernels: number of kernels used in the approximation (equals to the rank)
	%	- nSamples: number of random samples to be generated in randomization
	%	- tolerance: the maximum tolerable relative entropy between the original and approximated arrays
    %
    % Output:
	%	- randomizedInputDistribution [nTags * nStates]: the input probability distribution when the sum rate is maximized
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent randomized input combination probability distribution
	%	- randomizedRate [nTags * 1]: the individual tag rate when the sum rate is maximized
    %
    % Comment:
    %   - generate random samples of candidate tag probability vectors whose outer product expectation equals to the joint input distribution array
	%	- equivalent to a random search with guiduance on the correlation matrix of the distributions
	%	- the number of kernels and random samples are designable (performance-complexity tradeoff)
	%	- the candidates follows uniform distribution within a sphere bounded within the probability simplex
    %
    % Author & Date: Yang (i@snowztail.com), 22 Mar 3

	% * Declare default tolerance
	arguments
		jointArray;
		dmtc;
		nKernels = 4;
		nSamples = 5e2;
		tolerance = 1e-6;
	end

	% * Get data
	nStates = size(jointArray, 1);
	nTags = ndims_modified(jointArray);

	% * Initialization
	powerSet = power_set(1 : nTags);
	nCases = length(powerSet);

	% * Randomization
	isConverged = false;
	kernelVector = rand(nStates, nTags, nKernels);
	kernelVector = kernelVector ./ sum(kernelVector, 1);
	kernelArray = cell(nKernels, 1);
	while ~isConverged
		cvx_begin
			variable kernelCoefficient(nKernels, 1)
			for iKernel = 1 : nKernels
				kernelVectorCell = num2cell(kernelVector(:, :, iKernel), 1);
				kernelArray{iKernel} = kernelCoefficient(iKernel) * outer_product(kernelVectorCell{:});
			end
			referenceArray = sum(cat(nTags + 1, kernelArray{:}), nTags + 1);
			relativeEntropy = sum_nested(rel_entr(jointArray, referenceArray), 1 : nTags);
			minimize relativeEntropy
			subject to
				kernelCoefficient == simplex(nKernels);
		cvx_end

		for iTag = 1 : nTags
			cvx_begin
				variable objectVector(nStates, nKernels)
				for iKernel = 1 : nKernels
					kernelVectorCell = num2cell(kernelVector(:, :, iKernel), 1);
					kernelVectorCell{iTag} = objectVector(:, iKernel);
					kernelArray{iKernel} = kernelCoefficient(iKernel) * outer_product(kernelVectorCell{:});
				end
				referenceArray = sum(cat(nTags + 1, kernelArray{:}), nTags + 1);
				relativeEntropy = sum_nested(rel_entr(jointArray, referenceArray), 1 : nTags);
				minimize relativeEntropy
				subject to
					for iKernel = 1 : nKernels
						objectVector(:, iKernel) == simplex(nStates);
					end
			cvx_end
			kernelVector(:, iTag, :) = objectVector;
		end
		isConverged = abs(relativeEntropy) <= tolerance;
	end
	% * Generate random vectors with prescribed mean by kernel vectors
	randomizedRate = [];
	randomizedInputDistribution = [];
	for iKernel = 1 : nKernels
		for iSample = 1 : round(kernelCoefficient(iKernel) * nSamples)
			projectVector = zeros(nTags, nStates);
			for iTag = 1 : nTags
				% * Retrieve radius of the uniform random vector whose boundary lies within the probability simplex
				meanVector = kernelVector(:, iTag, iKernel);
				radiusBound(nStates + 1) = norm(meanVector) ^ 2;
				for iState = 1 : nStates
					radiusBound(iState) = abs(sum(meanVector(setdiff(1 : nStates, iState))) - 1) / sqrt(nStates - 1);
				end
				radius = min(radiusBound);
				% * Generate random vector uniformly distributed within the (nStates-dimensional) sphere
				isValid = false;
				while ~isValid
					randomVector = 2 * radius * rand(nStates, 1) + meanVector - radius;
					isValid = norm(randomVector - meanVector) <= radius;
				end
				projectVector(iTag, :) = transpose(randomVector + (1 - ones(1, nStates) * randomVector) / nStates * ones(nStates, 1));
			end
			% * Optimize rates with project vectors as input distributions
			rateBound = zeros(nCases, 1);
			cvx_begin
				variables rate(nTags, 1);
				expressions subSetRate(nCases, 1);
				for iCase = 1 : nCases
					subSetRate(iCase) = sum(rate(powerSet{iCase}));
					rateBound(iCase) = rate_bound_randomization(dmtc, powerSet{iCase}, projectVector);
				end
				% * Formulate problem
				maximize sum(rate)
				subject to
					for iCase = 1 : nCases
						0 <= subSetRate(iCase) <= rateBound(iCase);
					end
			cvx_end
			randomizedRate = cat(2, randomizedRate, rate);
			randomizedInputDistribution = cat(3, randomizedInputDistribution, projectVector);
		end
	end
	[~, randomizedIndex] = max(sum(randomizedRate, 1));
	randomizedRate = randomizedRate(:, randomizedIndex);
	randomizedInputDistribution = randomizedInputDistribution(:, :, randomizedIndex);
	randomizedEquivalentDistribution = prod(combination_distribution(randomizedInputDistribution), 1);
end
