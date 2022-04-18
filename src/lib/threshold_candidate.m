function [thresholdCandidate] = threshold_candidate(symbolRatio, equivalentChannel, noisePower, nBins, precoder, confidenceScore)
	% Function:
    %	- obtain finite decision threshold candidates (i.e., discretized output bin boundaries)
    %
    % Input:
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- noisePower: average noise power at the user
	%	- nBins: number of discretization bins over received signal
	%	- precoder [nTxs * 1]: transmit beamforming vector at the AP
    %	- confidenceScore: number of standard deviations from the mean (i.e., parameter k in Chebyshev's inequality)
    %
    % Output:
	%	- thresholdCandidate [1 * nBins + 1]: discretized output bin boundaries
    %
    % Comment:
	%	- obtain threshold candidates within empirical interval based on Chebyshev's inequality
	%	- increasing the number of candidates improves the backscatter achievable rate at the cost of (linear) computation complexity
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 18

	% * Declare default confidence score
	arguments
		symbolRatio;
		equivalentChannel;
		noisePower;
		nBins;
		precoder;
		confidenceScore = 10;
	end

	% * Compute the expected received power under all backscatter input combinations
	receivedPower = abs(equivalentChannel * precoder) .^ 2 + noisePower;

	% * Obtain threshold candidates within empirical interval based on Chebyshev's inequality
	lowerBound = symbolRatio * min(receivedPower) - confidenceScore * sqrt(symbolRatio * min(receivedPower));
	upperBound = symbolRatio * max(receivedPower) + confidenceScore * sqrt(symbolRatio * max(receivedPower));
	if lowerBound <= 0
		thresholdCandidate = [linspace(0, upperBound, nBins), inf];
	else
		thresholdCandidate = [0, linspace(lowerBound, upperBound, nBins - 1), inf];
	end
end
