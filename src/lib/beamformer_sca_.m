function [beamformer, dmtc, weightedSumRate] = beamformer_sca_(weight, symbolRatio, equivalentChannel, transmitPower, noisePower, equivalentDistribution, threshold, tolerance)
	% Function:
	%	- optimize the transmit beamformer to maximize the weighted sum primary-backscatter rate by successive convex approximation
    %
    % Input:
	%	- weight: the relative priority of the primary link
	%	- symbolRatio: the ratio of the backscatter symbol period over the primary symbol period
	%	- equivalentChannel [(nStates ^ nTags) * nTxs]: equivalent primary channel under each tag input combination
	%	- transmitPower: average transmit power at the AP
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
		transmitPower;
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
	beamformer = sqrt(transmitPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
% 	a = rand(1, nTxs) + 1i * rand(1, nTxs);
% 	beamformer = sqrt(transmitPower) * ctranspose(a) / norm(a);
	beamformerMatrix = beamformer * beamformer';

	% * Compute expected received power under all backscatter input combinations and update DMTC
	receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);
	dmtc = dmtc_local(symbolRatio, receivedPower, threshold);

	% * Initialize auxiliary variables associated with regularized incomplete gamma function
	epsilon = epsilon_auxiliary(receivedPower, threshold);
	xi = xi_auxiliary(symbolRatio, epsilon);
	[regularizedGamma, regularizedUpperGamma] = regularized_gamma(epsilon, xi);
	[xiLower, xiUpper] = deal(xi);
	[epsilonLower, epsilonUpper] = deal(epsilon);

	% * Iteratively update beamforming matrix by SCA
	isConverged = false;
	weightedSumRate = rate_weighted_sum(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformer, dmtc);
	while ~isConverged
		% * Update iteration index
		weightedSumRate_ = weightedSumRate;
		[receivedPower_, epsilon_, xi_, regularizedGamma_, regularizedUpperGamma_] = deal(receivedPower, epsilon, xi, regularizedGamma, regularizedUpperGamma);
		[xiLower_, xiUpper_] = deal(xi_);
		[epsilonLower_, epsilonUpper_] = deal(epsilon_);
% 		[xiLower_, xiUpper_] = deal(xiLower, xiUpper);
% 		[epsilonLower_, epsilonUpper_] = deal(epsilonLower, epsilonUpper);

		cvx_begin
			variable beamformerMatrix(nTxs, nTxs) hermitian semidefinite;
			variable regularizedGamma(nInputs, nOutputs) nonnegative;
			variables epsilonLower(nInputs, nOutputs + 1) epsilonUpper(nInputs, nOutputs + 1) xiLower(nInputs, nOutputs + 1) xiUpper(nInputs, nOutputs + 1);
			% variable epsilonLower(nInputs, nOutputs + 1) nonnegative;
			% variable epsilonUpper(nInputs, nOutputs + 1) nonnegative;
			% variable xiLower(nInputs, nOutputs + 1) nonnegative;
			% variable xiUpper(nInputs, nOutputs + 1) nonnegative;
			expressions primaryRate(nInputs, 1);
			expressions backscatterRateSca(nInputs, nOutputs);
			expressions z1(nInputs, nOutputs) z2(nInputs, nOutputs) c1(nInputs, nOutputs + 1) c2(nInputs, nOutputs + 1) c3(nInputs, nOutputs + 1) c4(nInputs, nOutputs + 1);
			% expressions regularizedGamma(nInputs, nOutputs);

			for iInput = 1 : nInputs
				primaryRate(iInput) = equivalentDistribution(iInput) * symbolRatio * log(1 + trace(equivalentChannel(iInput, :)' * equivalentChannel(iInput, :) * beamformerMatrix) / noisePower);
				for iOutput = 1 : nOutputs
					backscatterRateSca(iInput, iOutput) = equivalentDistribution(iInput) * regularizedGamma(iInput, iOutput) * log(regularizedGamma_(iInput, iOutput) / (equivalentDistribution * regularizedGamma_(:, iOutput)));
				end
			end
			receivedPower = received_power(equivalentChannel, noisePower, beamformerMatrix);

			maximize weight * sum(primaryRate) + (1 - weight) * sum(sum(backscatterRateSca))
			subject to
				beamformerMatrix == hermitian_semidefinite(nTxs);
				trace(beamformerMatrix) <= transmitPower;
				for iInput = 1 : nInputs
					for iOutput = 1 : nOutputs + 1
						% * epsilon
						c1(iInput, iOutput) = square(epsilonUpper(iInput, iOutput) - receivedPower(iInput)) - (epsilonUpper_(iInput, iOutput) ^ 2 + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * 2 * epsilonUpper_(iInput, iOutput)) ...
							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput));
						c2(iInput, iOutput) = square(epsilonLower(iInput, iOutput) + receivedPower(iInput)) - (epsilonLower_(iInput, iOutput) ^ 2 + (epsilonLower(iInput, iOutput) - epsilonLower_(iInput, iOutput)) * 2 * epsilonLower_(iInput, iOutput)) ...
							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput));
						c1(iInput, iOutput) <= - 2 * threshold(iOutput);
						c2(iInput, iOutput) <= 2 * threshold(iOutput);
% 						square(epsilonUpper(iInput, iOutput) - receivedPower(iInput)) - (epsilonUpper_(iInput, iOutput) ^ 2 + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * 2 * epsilonUpper_(iInput, iOutput)) ...
% 							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= - 2 * threshold(iOutput);
% 						square(epsilonLower(iInput, iOutput) + receivedPower(iInput)) - (epsilonLower_(iInput, iOutput) ^ 2 + (epsilonLower(iInput, iOutput) - epsilonLower_(iInput, iOutput)) * 2 * epsilonLower_(iInput, iOutput)) ...
% 							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= 2 * threshold(iOutput);
% 						square(epsilonLower(iInput, iOutput) - receivedPower(iInput)) - (epsilonLower_(iInput, iOutput) ^ 2 + (epsilonLower(iInput, iOutput) - epsilonLower_(iInput, iOutput)) * 2 * epsilonLower_(iInput, iOutput)) ...
% 							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= - 2 * threshold(iOutput);
% 						square(epsilonUpper(iInput, iOutput) + receivedPower(iInput)) - (epsilonUpper_(iInput, iOutput) ^ 2 + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * 2 * epsilonUpper_(iInput, iOutput)) ...
% 							- (receivedPower_(iInput) ^ 2 + (receivedPower(iInput) - receivedPower_(iInput)) * 2 * receivedPower_(iInput)) <= 2 * threshold(iOutput);
						% * xi
						c3(iInput, iOutput) = sum(pow_p(epsilonUpper(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
						c4(iInput, iOutput) = sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1)) + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 2) ./ factorial(0 : symbolRatio - 2));
						xiUpper(iInput, iOutput) >= c3(iInput, iOutput);
						xiLower(iInput, iOutput) <= c4(iInput, iOutput);
% 						xiUpper(iInput, iOutput) >= sum(pow_p(epsilonUpper(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
% 						xiLower(iInput, iOutput) <= sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1)) + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 2) ./ factorial(0 : symbolRatio - 2));
% 						xiLower(iInput, iOutput) >= sum(pow_p(epsilonUpper(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1));
% 						xiUpper(iInput, iOutput) <= sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 1) ./ factorial(0 : symbolRatio - 1)) + (epsilonUpper(iInput, iOutput) - epsilonUpper_(iInput, iOutput)) * sum(pow_p(epsilonUpper_(iInput, iOutput), 0 : symbolRatio - 2) ./ factorial(0 : symbolRatio - 2));
					end
				end

				for iInput = 1 : nInputs
					for iOutput = 1 : nOutputs
						% * regularized incomplete gamma
% 						z1(iInput, iOutput) = (exp(-epsilonUpper(iInput, iOutput) + (xiLower(iInput, iOutput) - xiLower_(iInput, iOutput)) / xiLower_(iInput, iOutput)) * xiLower_(iInput, iOutput)) ...
% 							- ((-epsilonLower(iInput, iOutput + 1) + log(xiUpper(iInput, iOutput + 1)) - log(regularizedUpperGamma_(iInput, iOutput + 1))) * regularizedUpperGamma_(iInput, iOutput + 1) + regularizedUpperGamma_(iInput, iOutput + 1));
						z2(iInput, iOutput) = ((-epsilonLower(iInput, iOutput) + log(xiUpper(iInput, iOutput)) - log(regularizedUpperGamma_(iInput, iOutput))) * regularizedUpperGamma_(iInput, iOutput) + regularizedUpperGamma_(iInput, iOutput)) ...
							- (exp(-epsilonUpper(iInput, iOutput + 1) + (xiLower(iInput, iOutput + 1) - xiLower_(iInput, iOutput + 1)) / xiLower_(iInput, iOutput + 1)) * xiLower_(iInput, iOutput + 1));
% 						regularizedGamma(iInput, iOutput) >= z1(iInput, iOutput);
						regularizedGamma(iInput, iOutput) <= z2(iInput, iOutput);
% 						regularizedGamma(iInput, iOutput) >= (exp(-epsilonUpper(iInput, iOutput) + (xiLower(iInput, iOutput) - xiLower_(iInput, iOutput)) / xiLower_(iInput, iOutput)) * xiLower_(iInput, iOutput)) ...
% 							- ((-epsilonLower(iInput, iOutput + 1) + log(xiUpper(iInput, iOutput + 1)) - log(regularizedUpperGamma_(iInput, iOutput + 1))) * regularizedUpperGamma_(iInput, iOutput + 1) + regularizedUpperGamma_(iInput, iOutput + 1));
% 						regularizedGamma(iInput, iOutput) <= ((-epsilonLower(iInput, iOutput) + log(xiUpper(iInput, iOutput)) - log(regularizedUpperGamma_(iInput, iOutput))) * regularizedUpperGamma_(iInput, iOutput) + regularizedUpperGamma_(iInput, iOutput)) ...
% 							- (exp(-epsilonUpper(iInput, iOutput + 1) + (xiLower(iInput, iOutput + 1) - xiLower_(iInput, iOutput + 1)) / xiLower_(iInput, iOutput + 1)) * xiLower_(iInput, iOutput + 1));
% 						regularizedGamma(iInput, iOutput) >= (exp(-epsilonLower(iInput, iOutput) + (xiUpper(iInput, iOutput) - xiUpper_(iInput, iOutput)) / xiUpper_(iInput, iOutput)) * xiUpper_(iInput, iOutput)) ...
% 							- ((-epsilonUpper(iInput, iOutput + 1) + log(xiLower(iInput, iOutput + 1)) - log(regularizedUpperGamma_(iInput, iOutput + 1))) * regularizedUpperGamma_(iInput, iOutput + 1) + regularizedUpperGamma_(iInput, iOutput + 1));
% 						regularizedGamma(iInput, iOutput) <= ((-epsilonUpper(iInput, iOutput) + log(xiLower(iInput, iOutput)) - log(regularizedUpperGamma_(iInput, iOutput))) * regularizedUpperGamma_(iInput, iOutput) + regularizedUpperGamma_(iInput, iOutput)) ...
% 							- (exp(-epsilonLower(iInput, iOutput + 1) + (xiUpper(iInput, iOutput + 1) - xiUpper_(iInput, iOutput + 1)) / xiUpper_(iInput, iOutput + 1)) * xiUpper_(iInput, iOutput + 1));
					end
				end

				for iInput = 1 : nInputs
					transpose(regularizedGamma(iInput, :)) == simplex(nOutputs);
% 					transpose(epsilonLower(iInput, :)) == nonnegative(nOutputs + 1);
% 					transpose(epsilonUpper(iInput, :)) == nonnegative(nOutputs + 1);
% 					transpose(xiLower(iInput, :)) == nonnegative(nOutputs + 1);
% 					transpose(xiUpper(iInput, :)) == nonnegative(nOutputs + 1);
				end
		cvx_end
		regularizedGamma2 = regularizedGamma;

		% * Update auxiliary variables
		epsilon = epsilon_auxiliary(receivedPower, threshold);
		xi = xi_auxiliary(symbolRatio, epsilon);
		[regularizedGamma, regularizedUpperGamma] = regularized_gamma(epsilon, xi);

		% * Update DMTC and compute actual weighted sum rate
		dmtc = dmtc_local(symbolRatio, receivedPower, threshold);
		[weightedSumRate, rate] = rate_weighted_sum_local(weight, symbolRatio, noisePower, equivalentDistribution, receivedPower, dmtc);

		% * Test convergence
		isConverged = abs(weightedSumRate - weightedSumRate_) <= tolerance;
	end

	% % ? NLOPT by fmincon
	% clearvars receivedPower regularizedGamma
	% for iInput = 1 : nInputs
	% 	% receivedPower(iInput) = @(W) (trace(equivalentChannel(iInput, :)' * equivalentChannel(iInput, :) * W) + noisePower);
	% 	conditionalEnergy = @(z, W) (z .^ (symbolRatio - 1) .* exp(-z ./ trace(equivalentChannel(iInput, :)' * equivalentChannel(iInput, :) * W) + noisePower)) ./ (trace(equivalentChannel(iInput, :)' * equivalentChannel(iInput, :) * W) + noisePower .^ symbolRatio .* gamma(symbolRatio));
	% 	for iOutput = 1 : nOutputs
	% 		regularizedGamma{iInput, iOutput} = @(W) integral(@(z) conditionalEnergy, threshold(iOutput), threshold(iOutput + 1));
	% 	end
	% end

	% backscatterRateNlopt = @(W) 0;
	% for iInput = 1 : nInputs
	% 	for iOutput = 1 : nOutputs
	% 		backscatterRateNlopt = @(W) backscatterRateNlopt(W) + equivalentDistribution(iInput) * regularizedGamma{iInput, iOutput}(W) * log(regularizedGamma{iInput, iOutput}(W) / (equivalentDistribution * regularizedGamma{:, iOutput}(W)));
	% 	end
	% end

	% fmincon(-backscatterRateNlopt, beamformerMatrix, )

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

function [weightedSumRate, rate] = rate_weighted_sum_local(weight, symbolRatio, noisePower, equivalentDistribution, receivedPower, dmtc)
	[nInputs, nOutputs] = deal(size(receivedPower, 1));
	% * Primary achievable rate
	primaryInformationFunction = symbolRatio * log(1 + receivedPower / noisePower);
	primaryRate = equivalentDistribution * primaryInformationFunction;

	% * Backscatter sum rate
	backscatterInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	backscatterRate = equivalentDistribution * sum(backscatterInformationFunction, 2);

	% * Weighted sum rate
	rate = [primaryRate; backscatterRate];
	weightedSumRate = [weight, 1 - weight] * rate;
end
