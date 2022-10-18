function [wsr, rate] = rate_weighted(snr, equivalentDistribution, dmac, weight)
	% Function:
    %	- compute weighted sum rate, primary rate, and total backscatter rate
    %
    % Input:
	%	- snr [nInputs x 1]: average receive signal-to-noise ratio per primary symbol for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- dmac [nInputs x nOutputs]: discrete memoryless thresholding multiple access channel whose input and output are tag state tuple
	%	- weight: relative priority of primary link
    %
    % Output:
	%	- wsr: weighted sum of primary rate and total backscatter rate
	%	- rate [2 x 1]: primary rate (nats/s/Hz) and total backscatter rate (nats/cu)
    %
    % Author & Date: Yang (i@snowztail.com), 22 Feb 22

	% * Set default weight
	arguments
		snr;
		equivalentDistribution;
		dmac;
		weight = 0;
	end

	primaryRate = equivalentDistribution' * log(1 + snr);
	backscatterRate = sum(entr(equivalentDistribution' * dmac) - equivalentDistribution' * entr(dmac));
	rate = [primaryRate; backscatterRate];
	wsr = weight * primaryRate + (1 - weight) * backscatterRate;
end
