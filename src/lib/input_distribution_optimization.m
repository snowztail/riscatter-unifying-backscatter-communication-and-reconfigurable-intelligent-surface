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
		tolerance = eps;
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
	subSetRate = cvx(zeros(nCases, 1));
	rateBound = cvx(zeros(nCases, 1));

	% * Optimize correlated input probability arrays
	cvx_begin
		% * Create variables and expressions
		variables rate(nTags, 1);
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
	weightedSumRate = sum(rate);

	% * Randomization
	nKernels = 4;
	isConverged = false;
	% referenceArray = 0;
	kernelVector = rand(nStates, nTags, nKernels);
	kernelVector = kernelVector ./ sum(kernelVector, 1);
	kernelArray = cell(nKernels, 1);
	while ~isConverged
		cvx_begin
			variable kernelCoefficient(nKernels, 1)
			for iKernel = 1 : nKernels
				pCell = num2cell(kernelVector(:, :, iKernel), 1);
				% referenceArray = referenceArray + kernelCoefficient(iKernel) * outer_product(pCell{:});
				kernelArray{iKernel} = kernelCoefficient(iKernel) * outer_product(pCell{:});
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
					pCell = num2cell(kernelVector(:, :, iKernel), 1);
					pCell{iTag} = objectVector(:, iKernel);
					kernelArray{iKernel} = kernelCoefficient(iKernel) * outer_product(pCell{:});
				end
				referenceArray = sum(cat(nTags + 1, kernelArray{:}), nTags + 1);
				relativeEntropy = sum_nested(rel_entr(jointArray, referenceArray), 1 : nTags);
				minimize relativeEntropy
				subject to
					for iKernel = 1 : nKernels
						objectVector(:, iKernel) == simplex(nStates);
					end
			cvx_end
% 			kernelVector(:, iTag, :) = permute(objectVector, [1 3 2]);
			kernelVector(:, iTag, :) = objectVector;
		end
		isConverged = abs(relativeEntropy) <= eps;
	end
end
