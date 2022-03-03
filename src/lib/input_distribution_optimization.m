function [weightedSumRate] = input_distribution_optimization(nTags, dmtc, weight, symbolRatio, snr, tolerance)
	% Function:
	%	- optimize the tag input distribution to characterize the capacity region of user and tags
    %
    % Input:
	%	- nTags: number of tags
    %   - dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weight [2 * 1]: the relative priority of the primary and backscatter links
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- snr [(nStates ^ nTags) * 1]: signal-to-noise ratio of the primary link corresponding to to each input letter combination
    %   - tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- inputDistribution [nTags * nStates]: input probability distribution
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (bps/Hz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (bpcu)
    %
    % Comment:
    %   - Obtain the joint input probability matrix by optimization (corresponding to the optimal input with full transmit cooperation)
	%	- Extract rank-1 solution from the joint input probability matrix by randomization or marginalization
	%	- the discrete memoryless MAC is given in joint (equivalent point-to-point) form P(y | x_1, ..., x_K), instead of marginal form p(y | x_k)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 23

	% * Declare default tolerance
	arguments
		nTags;
		dmtc;
		weight;
		symbolRatio;
		snr;
		tolerance = 1e-6;
	end

	% * Ensure non-zero channel transition probability
	dmtc(dmtc < eps) = eps;
	dmtc = dmtc ./ sum(dmtc, 2);

	% * Get data
	nStates = nthroot(size(dmtc, 1), nTags);

	% * Initialization
	[powerSet, complementSet] = power_set(1 : nTags);
	nCases = length(powerSet);
	inputArraySet = cell(nCases, 1);

	% * Optimize correlated input probability arrays
	cvx_begin
		% * Create variables and expressions
		variables rate(nTags, 1);
		expressions subSetRate(nCases, 1) rateBound(nCases, 1);
		for iCase = 1 : nCases
			inputArraySet{iCase} = strcat('inputArray', sprintf('%d', powerSet{iCase}));
			eval(['variable', ' ', inputArraySet{iCase}, '(', num2str(repmat(nStates,[1, length(powerSet{iCase})])), ')']);
			subSetRate(iCase) = sum(rate(powerSet{iCase}));
		end

		% * Obtain upper bounds of sum rate on all subsets
		jointArray = eval(inputArraySet{end});
		for iCase = 1 : nCases
			if iCase < nCases
				complementArray = eval(inputArraySet{nCases - iCase});
			else
				complementArray = [];
			end
			rateBound(iCase) = rate_bound(dmtc, powerSet{iCase}, jointArray, complementArray);
		end

		% * Formulate problem
		maximize sum(rate)
		subject to
			for iCase = 1 : nCases
				0 <= subSetRate(iCase) <= rateBound(iCase);
				permute(sum_nested(eval(inputArraySet{end}), complementSet{iCase}), [powerSet{iCase}, complementSet{iCase}]) == eval(inputArraySet{iCase});
			end
			eval(inputArraySet{end}) == semidefinite(size(eval(inputArraySet{end})));
			sum_nested(eval(inputArraySet{end}), powerSet{end}) == 1;
	cvx_end
	relaxedRate = rate;
	weightedSumRate = sum(relaxedRate);

	% * Randomization
	nKernels = 4;
	nSamples = 10;
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
				expressions subSetRate(nCases, 1) rateBound(nCases, 1);
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
end
