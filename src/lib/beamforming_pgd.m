function [beamforming] = beamforming_pgd(symbolRatio, weight, transmitPower, noisePower, equivalentChannel, backwardChannel, equivalentDistribution, threshold, tolerance, alpha, beta)
	% Function:
	%	- optimize beamforming vector by projected gradient descent with step size determined by backtracking line search
    %
    % Input:
	%	- symbolRatio: backscatter/primary symbol duration ratio
	%	- weight: relative priority of primary link
	%	- transmitPower: average transmit power
	%	- noisePower: average noise power
	%	- equivalentChannel [nTxs x nInputs]: equivalent primary channel for each tag state tuple
	%	- equivalentDistribution [nInputs x 1]: equivalent single-source distribution for each tag input distribution tuple
	%	- threshold [1 x (nOutputs + 1)]: boundaries of quantization bins or decision regions (including 0 and Inf)
    %	- tolerance: minimum rate gain ratio per iteration
	%	- alpha: fraction of acceptable decrease predicted by linear extrapolation
	%	- beta: ratio of step size reduction in backtracking line search
    %
    % Output:
	%	- beamforming [nTxs x 1]: transmit beamforming vector
    %
	% Comment:
	%	- complex gradient w.r.t. beamforming vector is obtained in closed form
	%	- use backtracking line search to accelerate convergence
	%	- when primary link is prioritized, MRT to (ergodic) equivalent channel is nearly optimal; this may not hold when backscatter links are prioritized
	%
    % Author & Date: Yang (i@snowztail.com), 22 May 16

	% * Set default tolerance and parameters of backtracking line search
	arguments
		symbolRatio;
		weight;
		transmitPower;
		noisePower;
		equivalentChannel;
		backwardChannel;
		equivalentDistribution;
		threshold;
		tolerance = 1e-3;
		alpha = 1e-2;
		beta = 0.5;
	end

	% ! Initialize beamforming by previous solution
	persistent initializer

	% * No previous solution, initialize randomly
% 	if isempty(initializer)
	if true
% 		ric = weight * equivalentChannel * equivalentDistribution + (1 - weight) * sum(cascadedChannel, 2);
		ric = weight * equivalentChannel * equivalentDistribution + (1 - weight) * backwardChannel * equivalentDistribution;
% 		ric = rand(size(equivalentChannel, 1), 1) + 1i * rand(size(equivalentChannel, 1), 1);
		initializer.beamforming = sqrt(transmitPower) * ric / norm(ric);
	end

	% * Apply initializer
	beamforming = initializer.beamforming;

	% * Construct DMTC and recover i/o mapping
	receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
	snr = receivePower / noisePower;
	[~, sortIndex] = sort(receivePower);
	dmac = dmc_integration(symbolRatio, receivePower, threshold);
	dmac(:, sortIndex) = dmac;
	wsr = rate_weighted(weight, snr, equivalentDistribution, dmac);

	% * Projected gradient descent
	isConverged = false;
	while ~isConverged
		% * Compute local gradient
		gradient = gradient_local(symbolRatio, weight, noisePower, equivalentChannel, equivalentDistribution, threshold, beamforming);

		% * Optimize step size by backtracking line search
		step = 1;
		beamformingPgd = beamforming_projection(transmitPower, beamforming + step * gradient);
		receivePowerPgd = abs(equivalentChannel' * beamformingPgd) .^ 2 + noisePower;
		snrPgd = receivePowerPgd / noisePower;
		dmacPgd = dmc_integration(symbolRatio, receivePowerPgd, threshold);
		wsrPgd = rate_weighted(weight, snrPgd, equivalentDistribution, dmacPgd);
		while wsr > wsrPgd + alpha * step * norm(gradient) ^ 2 && step > eps
			step = beta * step;
			beamformingPgd = beamforming_projection(transmitPower, beamforming + step * gradient);
			receivePowerPgd = abs(equivalentChannel' * beamformingPgd) .^ 2 + noisePower;
			snrPgd = receivePowerPgd / noisePower;
			dmacPgd = dmc_integration(symbolRatio, receivePowerPgd, threshold);
			wsrPgd = rate_weighted(weight, snrPgd, equivalentDistribution, dmacPgd);
		end

		% * Test convergence (gradient can be non-zero due to norm constraint)
		isConverged = (wsrPgd - wsr) / wsr <= tolerance || any(dmacPgd < tolerance, 'all') || isnan(wsrPgd);
		beamforming = beamformingPgd;
		wsr = wsrPgd;
	end

	% * Update initializer
	initializer.beamforming = beamforming;
end


function [gradient] = gradient_local(symbolRatio, weight, noisePower, equivalentChannel, equivalentDistribution, threshold, beamforming)
	% * Get data
	[nInputs, nOutputs] = deal(size(equivalentChannel, 2));
	nTxs = size(beamforming, 1);

	% * Evaluate expressions
	receivePower = abs(equivalentChannel' * beamforming) .^ 2 + noisePower;
	dmac = dmc_integration(symbolRatio, receivePower, threshold);

	% * Compute derivative and gradient (loops can be further optimized)
	dDmac = zeros(nTxs, nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			coefficient = ((symbolRatio - 1 : -1 : 1) - threshold(iOutput + 1) / receivePower(iInput)) ./ factorial(symbolRatio - 1 : -1 : 1);
			dDmac(:, iInput, iOutput) = equivalentChannel(:, iInput) * equivalentChannel(:, iInput)' * beamforming / receivePower(iInput) ...
				* (threshold(iOutput + 1) * exp(-threshold(iOutput + 1) / receivePower(iInput)) * (polyval(coefficient, threshold(iOutput + 1) / receivePower(iInput)) - 1) ...
					- threshold(iOutput) * exp(-threshold(iOutput) / receivePower(iInput)) * (polyval(coefficient, threshold(iOutput) / receivePower(iInput)) - 1));
		end
	end
	dWsr = zeros(nTxs, 1);
	for iInput = 1 : nInputs
		dWsr = dWsr + weight * equivalentDistribution(iInput) * equivalentChannel(:, iInput) * equivalentChannel(:, iInput)' * beamforming / receivePower(iInput);
		for iOutput = 1 : nOutputs
			dWsr = dWsr + (1 - weight) * equivalentDistribution(iInput) * ((log(dmac(iInput, iOutput) / (equivalentDistribution' * dmac(:, iOutput))) + 1) * dDmac(:, iInput, iOutput) ...
				- dmac(iInput, iOutput) * dDmac(:, :, iOutput) * equivalentDistribution / (equivalentDistribution' * dmac(:, iOutput)));
		end
	end
	gradient = 2 * dWsr;
end

function [beamforming] = beamforming_projection(transmitPower, beamforming)
	beamforming = sqrt(transmitPower) * beamforming / max(sqrt(transmitPower), norm(beamforming));
end
