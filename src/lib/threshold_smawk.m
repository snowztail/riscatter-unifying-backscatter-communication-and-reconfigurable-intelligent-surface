function [threshold] = threshold_smawk(equivalentDistribution, thresholdDomain, binDmc)
	% Function:
	%	- group received energy bins into convex adjoint decision regions by dynamic programming accelerated by SMAWK algorithm
    %
    % Input:
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- thresholdDomain [1 x (nBins + 1)]: boundaries of quantized energy bins as domain of decision thresholds	%	- binDmc [nInputs x nBins]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins
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

	% ! Pre-compute quantization cost to accelerate runtime
	cost = zeros(nBins);
	for iBin = 1 : nBins
		for jBin = iBin : nBins
			cost(iBin, jBin) = cost_quantization(iBin : jBin, jointDistribution, outputDistribution);
		end
	end

	% * Initialize cost functions
	dp = zeros(nBins, nOutputs);
	sol = zeros(nBins, nOutputs);
	for iBin = 1 : nBins
		% dp(iBin, 1) = cost_quantization(1 : iBin, jointDistribution, outputDistribution);
		dp(iBin, 1) = cost(1, iBin);
	end

	% * Generate monotone matrix
	for iOutput = 2 : nOutputs
		% d = Inf(nBins - nOutputs + 1);
		% for iRow = 1 : nBins - nOutputs + 1
		% 	for iColumn = 1 : iRow
		% 		% d(iRow, iColumn) = dp(iColumn - 2 + iOutput, iOutput - 1) + cost_quantization(iColumn - 1 + iOutput : iRow - 1 + iOutput, jointDistribution, outputDistribution);
		% 		d(iRow, iColumn) = dp(iColumn - 2 + iOutput, iOutput - 1) + cost(iColumn - 1 + iOutput, iRow - 1 + iOutput);
		% 	end
		% end

		d = Inf(nBins - nOutputs + 1);
		for iRow = 1 : nBins - nOutputs + 1
			% d(iRow, 1 : iRow) = dp(1 - 2 + iOutput : iRow - 2 + iOutput, iOutput - 1) + cost(1 - 1 + iOutput : iRow - 1 + iOutput, iRow - 1 + iOutput);
			d(iRow, 1 : iRow) = dp(-1 + iOutput : iRow - 2 + iOutput, iOutput - 1) + cost(iOutput : iRow - 1 + iOutput, iRow - 1 + iOutput);
		end

		% * Retrieve leftmost minima position by SMAWK
		[p] = smawk(d, 1 : nBins - nOutputs + 1, 1 : nBins - nOutputs + 1);

		% * Get sol and dp
		for iBin = nBins - nOutputs + iOutput : -1 : iOutput
			sol(iBin, iOutput) = p(iBin - iOutput + 1) - 2 + iOutput;
			% dp(iBin, iOutput) = dp(sol(iBin, iOutput), iOutput - 1) + cost_quantization(sol(iBin, iOutput) + 1 : iBin, jointDistribution, outputDistribution);
			dp(iBin, iOutput) = dp(sol(iBin, iOutput), iOutput - 1) + cost(sol(iBin, iOutput) + 1, iBin);
		end
	end

	% * Recursively generate thresholds
	index(nOutputs + 1, 1) = nBins;
	for iOutput = nOutputs : -1 : 1
		index(iOutput) = sol(index(iOutput + 1), iOutput);
	end
	threshold = thresholdDomain(index + 1);
end


function [p] = smawk(d, r, c)
	p = zeros(1, length(r));
	[c] = reduce(d, r, c);
	if length(r) == 1
		p = c;
	else
		[p(2 : 2 : end)] = smawk(d, r(2 : 2 : end), c);
		j = 1;
		for i = 1 : 2 : length(r)
			p(i) = c(j);
			if i < length(r)
				u = p(i + 1);
			else
				u = Inf;
			end
			while j <= length(r) && c(j) <= u
				if d(r(i), c(j)) < d(r(i), p(i))
					p(i) = c(j);
				end
				j = j + 1;
			end
			j = j - 1;
		end
	end
end

function [c] = reduce(d, r, c)
	i = 1;
	while length(r) < length(c)
		if d(r(i), c(i)) <= d(r(i), c(i + 1))
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
