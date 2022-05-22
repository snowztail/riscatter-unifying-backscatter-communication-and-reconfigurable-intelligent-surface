function [threshold] = threshold_smawk(equivalentDistribution, quantizationLevel, binDmc)
	% Function:
	%	- group received energy bins into convex adjoint decision regions by dynamic programming accelerated by SMAWK algorithm
    %
    % Input:
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- quantizationLevel [1 x (nBins + 1)]: boundaries of quantized energy bins
	%	- binDmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
    %
    % Output:
	%	- threshold [1 x (nOutputs + 1)]: boundaries of decision regions (including 0 and Inf)
    %
    % Comment:
	%	- quantization cost function satisfies Quadrangle Inequality (QI) that enables SMAWK algorithm
	%	- assume input and output alphabe size are same
    %
	% Reference:
	%	- X. He, K. Cai, W. Song, and Z. Mei, “Dynamic programming for sequential deterministic quantization of discrete memoryless channels,” IEEE Transactions on Communications, vol. 69, no. 6, pp. 3638–3651, 2021.
	%
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09


	% * Get data
	[nInputs, nBins] = size(binDmc);
	nOutputs = nInputs;

	% * Obtain relevant distributions
	jointDistribution = equivalentDistribution .* binDmc;
	outputDistribution = equivalentDistribution' * binDmc;

	% * Initialize cost functions
	dp = zeros(nBins, nOutputs);
	sol = zeros(nBins, nOutputs);
	for iBin = 1 : nBins
		dp(iBin, 1) = cost_quantization(1 : iBin, jointDistribution, outputDistribution);
	end

	% * Generate monotone matrix
	for iOutput = 2 : nOutputs
		D = zeros(nBins - nOutputs + 1);
		for iRow = 1 : nBins - nOutputs + 1
			for iColumn = 1 : nBins - nOutputs + 1
				if iRow >= iColumn
					D(iRow, iColumn) = dp(iColumn - 2 + iOutput, iOutput - 1) + cost_quantization(iColumn - 1 + iOutput : iRow - 1 + iOutput, jointDistribution, outputDistribution);
				else
					D(iRow, iColumn) = Inf;
				end
			end
		end

		% * Retrieve leftmost minima position by SMAWK
		[r, c] = deal(1 : nBins - nOutputs + 1);
		[p] = smawk(D, r, c);

		% * Get sol and dp
		for iBin = nBins - nOutputs + iOutput : -1 : iOutput
			sol(iBin, iOutput) = p(iBin - iOutput + 1) - 2 + iOutput;
			dp(iBin, iOutput) = dp(sol(iBin, iOutput), iOutput - 1) + cost_quantization(sol(iBin, iOutput) + 1 : iBin, jointDistribution, outputDistribution);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : -1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = quantizationLevel(index + 1);
end


function [p] = smawk(matrix, r, c)
	p = zeros(1, length(r));
	[c] = reduce(matrix, r, c);
	if length(r) == 1
		p = c;
	else
		[p(2 : 2 : end)] = smawk(matrix, r(2 : 2 : end), c);
		j = 1;
		for i = 1 : 2 : length(r)
			p(i) = c(j);
			if i < length(r)
				u = p(i + 1);
			else
				u = Inf;
			end
			while j <= length(r) && c(j) <= u
				if matrix(r(i), c(j)) < matrix(r(i), p(i))
					p(i) = c(j);
				end
				j = j + 1;
			end
			j = j - 1;
		end
	end
end

function [c] = reduce(matrix, r, c)
	i = 1;
	while length(r) < length(c)
		if matrix(r(i), c(i)) <= matrix(r(i), c(i + 1))
			if i < length(r)
				i = i + 1;
			elseif i == length(r)
				c(i + 1) = [];
			end
		else
			c(i) = [];
			if i > 1
				i = i - 1;
			end
		end
	end
end
