function [beamformer, dmtc, weightedSumRate] = beamformer_newton(weight, symbolRatio, equivalentChannel, transmitPower, noisePower, equivalentDistribution, threshold, tolerance)


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

	% * Initialize beamformer
	beamformer = sqrt(transmitPower) * ctranspose(equivalentDistribution * equivalentChannel) / norm(equivalentDistribution * equivalentChannel);
% 	beamformer = conj(beamformer);
	beamformer = rand(size(beamformer)) + 1i * rand(size(beamformer));
	beamformerComponent_ = zeros(nTxs, 2);
	beamformerComponent_(:, 1) = real(beamformer);
	beamformerComponent_(:, 2) = imag(beamformer);
% 	[receivedPower] = received_power(equivalentChannel, noisePower, beamformerComponent);

	options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
	f = @(beamformerComponent) -weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformerComponent_, threshold);
% 	[beamformerComponent, wsr] = fminunc(f, beamformerComponent_, options);
	[beamformerComponent, wsr,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(f, beamformerComponent_, [], [], [], [], [], [], [], options);
end


function [receivedPower] = received_power(equivalentChannel, noisePower, beamformerComponent)
	nInputs = size(equivalentChannel, 1);
	receivedPower = zeros(nInputs, 1);
	for iInput = 1 : nInputs
		channelComponent = zeros(size(transpose(beamformerComponent)));
		channelComponent(1, :) = real(equivalentChannel(iInput, :));
		channelComponent(2, :) = - imag(equivalentChannel(iInput, :));
% 		receivedPower(iInput) = channelComponent(1, :) * beamformerComponent(:, 1) + channelComponent(2, :) * beamformerComponent(:, 2) + 1i * (channelComponent(1, :) * beamformerComponent(:, 2) - channelComponent(2, :) * beamformerComponent(:, 1));
		receivedPower(iInput) = (channelComponent(1, :) * beamformerComponent(:, 1) + channelComponent(2, :) * beamformerComponent(:, 2)) ^ 2 + (channelComponent(1, :) * beamformerComponent(:, 2) - channelComponent(2, :) * beamformerComponent(:, 1)) ^ 2 + noisePower;
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


function [weightedSumRate] = weighted_sum_rate_local(weight, symbolRatio, equivalentChannel, noisePower, equivalentDistribution, beamformerComponent, threshold)
	[nInputs, nOutputs] = deal(size(equivalentDistribution, 2));
	receivedPower = received_power(equivalentChannel, noisePower, beamformerComponent);
	dmtc = dmtc_local(symbolRatio, receivedPower, threshold);
	primaryRate = equivalentDistribution * symbolRatio * log(1 + receivedPower / noisePower);
	backscatterInformationFunction = zeros(nInputs, nOutputs);
	for iInput = 1 : nInputs
		for iOutput = 1 : nOutputs
			backscatterInformationFunction(iInput, iOutput) = dmtc(iInput, iOutput) * log(dmtc(iInput, iOutput) / (equivalentDistribution * dmtc(:, iOutput)));
		end
	end
	backscatterRate = equivalentDistribution * sum(backscatterInformationFunction, 2);
	weightedSumRate = [weight, 1 - weight] * [primaryRate; backscatterRate];
end
