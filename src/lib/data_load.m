Result(nVariables, nWeights, nInstances) = struct('weight', [], 'rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

instanceSet = 1 : nInstances;
for iInstance = 1 : nInstances
	try
		Result(:, :, iInstance) = load(strcat(directory, 'instance_', num2str(iInstance)), 'Result').Result;
	catch
		instanceSet(instanceSet == iInstance) = [];
	end
end
Result = Result(:, :, instanceSet);
