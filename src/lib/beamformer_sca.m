function [beamformer, dmtc, weightedSumRate] = beamformer_sca(weight, symbolRatio, equivalentChannel, txPower, noisePower, equivalentDistribution, threshold, tolerance)
	% Function:
	%	- optimize the transmit beamformer to maximize the weighted sum primary-backscatter rate by successive convex approximation
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent AP-user channels under all backscatter input combinations
	%	- txPower: average transmit power at the AP
	%	- noisePower: average noise power at the user
	%	- equivalentDistribution [1 * (nStates ^ nTags)]: equivalent input combination probability distribution
	%	- threshold [1 * (nOutputs + 1)]: boundaries of decision regions
    %	- tolerance: minimum rate gain per iteration
    %
    % Output:
	%	- beamformer [nTxs * 1]: transmit beamforming vector at the AP
    %	- dmtc [(nStates ^ nTags) * nOutputs]: the transition probability matrix of the backscatter discrete memoryless thresholding MAC
	%	- weightedSumRate: weighted sum of primary rate and total backscatter rate
	%		- primaryRate: the achievable rate for the primary link (nats per second per Hertz)
	%		- backscatterRate: the achievable sum rate for the backscatter link (nats per channel use)
    %
    % Comment:
    %	- iteratively update the beamforming matrix until convergence
	%	- extract the best beamformer by Gaussian randomization
    %
    % Author & Date: Yang (i@snowztail.com), 22 Apr 19

	% * Declare default tolerance
	arguments
		weight;
		symbolRatio;
		equivalentChannel;
		txPower;
		noisePower;
		equivalentDistribution;
		threshold;
		tolerance = 1e-6;
	end

	% ! Use finite threshold to evaluate incomplete gamma function by series representation
	threshold(end) = 10 * threshold(end - 1);

	% * Get data
	[nInputs, nOutputs] = deal(size(equivalentChannel, 1));
	nTxs = size(equivalentChannel, 2);

	% * Initialize beamformer matrix
	beamformer = sqrt(txPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
	beamformerMatrix = beamformer * beamformer';

	% * Compute expected received power under all backscatter input combinations and update DMTC
	receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);
	dmtc = dmtc_local(symbolRatio, receivedPower, threshold);

	% ? Sort the expected received power and arrange the equivalent distribution correspondingly
	% receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);
	% [sorted, sortIndex] = sort(receivedPower);
	% sortedDistribution = equivalentDistribution(sortIndex);

	% * Initialize auxiliary variables associated with regularized incomplete gamma function
	epsilon = epsilon_auxiliary(receivedPower, threshold);
	xi = xi_auxiliary(symbolRatio, epsilon);
	[regularizedGamma, regularizedUpperGamma] = regularized_gamma(epsilon, xi);

	% * Iteratively update beamforming matrix by SCA
	isConverged = false;
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
	while ~isConverged
		% * Update iteration index
		weightedSumRate_ = weightedSumRate;
		[receivedPower_, epsilon_, xi_, regularizedGamma_, regularizedUpperGamma_] = deal(receivedPower, epsilon, xi, regularizedGamma, regularizedUpperGamma);
		[xiLower_, xiUpper_] = deal(xi_);
		[epsilonLower_, epsilonUpper_] = deal(epsilon_);

		cvx_begin
		cvx_solver mosek
			variable beamformerMatrix(nTxs, nTxs) hermitian semidefinite;
			variables epsilonLower(nInputs, nOutputs + 1) epsilonUpper(nInputs, nOutputs + 1);
			variables xiLower(nInputs, nOutputs + 1) xiUpper(nInputs, nOutputs + 1);
			variable regularizedGamma(nInputs, nOutputs) nonnegative;
			expressions backscatterRateSca;

			for iInput = 1 : nInputs
				for iOutput = 1 : nOutputs
					% * Marginal entropy
					backscatterRateSca = backscatterRateSca + equivalentDistribution(iInput) * (regularizedGamma_(iInput, iOutput) * log(regularizedGamma_(iInput, iOutput)) + (regularizedGamma(iInput, iOutput) - regularizedGamma_(iInput, iOutput)) * (log(regularizedGamma_(iInput, iOutput)) + 1));
				end
				% * Conditional entropy
				backscatterRateSca = backscatterRateSca + entr(equivalentDistribution * regularizedGamma(:, iOutput));
			end
			receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);

			maximize backscatterRateSca
			subject to
				for iInput = 1 : nInputs
					for iOutput = 1 : nOutputs + 1
						% * epsilon
						square(epsilonUpper(iInput, iOutput) - receivedPower(iInput)) - (epsilonUpper_(iInput, iOutput) ^ 2 + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * 2 * epsilonUpper_(iInput, iOutput)) ...
							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= - 2 * threshold(iOutput);
						square(epsilonLower(iInput, iOutput) + receivedPower(iInput)) - (epsilonLower_(iInput, iOutput) ^ 2 + (epsilonLower(iInput, iOutput) - epsilonLower_(iInput, iOutput)) * 2 * epsilonLower_(iInput, iOutput)) ...
							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= 2 * threshold(iOutput);
						% * xi
						xiUpper(iInput, iOutput) >= sum(pow_p(epsilonUpper(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
						xiLower(iInput, iOutput) <= sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1)) + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 2) ./ factorial(0 : symbolRatio - 2));
					end
				end

				for iInput = 1 : nInputs
					for iOutput = 1 : nOutputs
						% * regularized incomplete gamma
						regularizedGamma(iInput, iOutput) >= (exp(-epsilonUpper(iInput, iOutput) + (xiLower(iInput, iOutput) - xiLower_(iInput, iOutput)) / xiLower_(iInput, iOutput)) * xiLower_(iInput, iOutput)) ...
							- ((-epsilonLower(iInput, iOutput + 1) + log(xiUpper(iInput, iOutput + 1)) - log(regularizedUpperGamma_(iInput, iOutput + 1))) * regularizedUpperGamma_(iInput, iOutput + 1) + regularizedUpperGamma_(iInput, iOutput + 1));
						regularizedGamma(iInput, iOutput) <= ((-epsilonLower(iInput, iOutput) + log(xiUpper(iInput, iOutput)) - log(regularizedUpperGamma_(iInput, iOutput))) * regularizedUpperGamma_(iInput, iOutput) + regularizedUpperGamma_(iInput, iOutput)) ...
							- (exp(-epsilonUpper(iInput, iOutput + 1) + (xiLower(iInput, iOutput + 1) - xiLower_(iInput, iOutput + 1)) / xiLower_(iInput, iOutput + 1)) * xiLower_(iInput, iOutput + 1));
					end
				end
		cvx_end

		% * Compute actual weighted sum rate
		weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end

end


function [dmtc] = dmtc_local(symbolRatio, receivedPower, threshold)
	nInputs = size(receivedPower, 1);
	nOutputs = size(threshold, 2) - 1;
	dmtc = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		conditionalEnergy = @(z) (z .^ (symbolRatio - 1) .* exp(-z ./ receivedPower(iInput))) ./ (receivedPower(iInput) .^ symbolRatio .* gamma(symbolRatio));
		for iOutput = 1 : nOutputs
			dmtc(iInput, iOutput) = integral(conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
		end
	end
end

function [receivedPower] = received_power(equivalentChannel, noisePower, beamformerMatrix)
	nInputs = size(equivalentChannel, 1);
	if isa(beamformerMatrix, 'cvx')
		receivedPower = cvx(zeros(nInputs, 1));
	else
		receivedPower = zeros(nInputs, 1);
	end
	for iInput = 1 : nInputs
		receivedPower(iInput) = real(trace((equivalentChannel(iInput, :)' * equivalentChannel(iInput, :)) * beamformerMatrix)) + noisePower;
	end
end

function [epsilon] = epsilon_auxiliary(receivedPower, threshold)
	nInputs = size(receivedPower, 1);
	nOutputs = size(threshold, 2) - 1;
	epsilon = zeros(nInputs, nOutputs + 1);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs + 1
			epsilon(iInput, iOutput) = threshold(iOutput) / receivedPower(iInput);
		end
	end
end

function [xi] = xi_auxiliary(symbolRatio, epsilon)
	nInputs = size(epsilon, 1);
	nOutputs = size(epsilon, 2) - 1;
	xi = zeros(nInputs, nOutputs + 1);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs + 1
			xi(iInput, iOutput) = sum(epsilon(iInput, iOutput) .^ (0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
		end
	end
end

function [regularizedGamma, regularizedUpperGamma] = regularized_gamma(epsilon, xi)
	nInputs = size(epsilon, 1);
	nOutputs = size(epsilon, 2) - 1;
	regularizedUpperGamma = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs + 1
			regularizedUpperGamma(iInput, iOutput) = exp(-epsilon(iInput, iOutput)) * xi(iInput, iOutput);
		end
	end
	regularizedGamma = - diff(regularizedUpperGamma, 1, 2);
end
