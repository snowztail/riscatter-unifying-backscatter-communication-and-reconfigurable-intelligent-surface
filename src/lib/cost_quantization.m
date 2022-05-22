function [quantizationCost] = cost_quantization(binIndex, jointDistribution, outputDistribution)
	% Function:
	%	- compute quantization cost of grouping indexed energy bins into a decision region
    %
    % Input:
	%	- jointDistribution [nInputs x nBins]: joint probability distribution of input tag distribution tuples and output energy bins
	%	- outputDistribution [1 x nBins]: probability distribution of output energy bins
    %
    % Output:
	%	- quantizationCost: cost of grouping indexed energy bins into a decision region
    %
    % Comment:
	%	- quantization cost evaluates loss of mutual information
	%
    % Author & Date: Yang (i@snowztail.com), 22 May 21

	backwardDistribution = sum(jointDistribution(:, binIndex), 2) / sum(outputDistribution(binIndex), 2);
	quantizationCost = sum(entr(backwardDistribution)) * sum(outputDistribution(binIndex));
end
