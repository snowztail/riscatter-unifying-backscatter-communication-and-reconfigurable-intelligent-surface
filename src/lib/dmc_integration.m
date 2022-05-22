function [dmc] = dmc_integration(symbolRatio, receivePower, threshold)
	% Function:
	%	- construct discrete memoryless channel for energy detection based on conditional energy p.d.f. and decision thresholds
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- receivePower [nInputs x 1]: average receive power per primary symbol for each tag state tuple
	%	- threshold [1 x (nOutputs + 1)]: boundaries of quantization bins or decision regions (including 0 and Inf)
    %
    % Output:
	%	- dmc [nInputs x nOutputs]: discrete memoryless channel whose input is tag state tuple and output is (high-resolution) quantized energy bins or tag state tuple
    %
    % Comment:
    %	- for each tag state tuple, received power per primary symbol follows Erlang distribution
	%	- DMC depends on bin boundaries or decision thresholds
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 09

	% * Get data
	nInputs = size(receivePower, 1);
	nOutputs = size(threshold, 2) - 1;
	% [sortPower, sortIndex] = sort(receivePower);

	% * Construct DMC based on conditional energy p.d.f. and decision thresholds
	dmc = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		channelDistribution = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivePower(iInput))) ./ (receivePower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			dmc(iInput, iOutput) = integral(channelDistribution, threshold(iOutput), threshold(iOutput + 1));
		end
	end

	% % * Desort to recover i/o mapping
	% dmc(:, sortIndex) = dmc;
end
