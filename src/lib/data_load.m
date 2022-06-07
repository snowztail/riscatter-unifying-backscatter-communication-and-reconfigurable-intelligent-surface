Result(nVariables, nInstances) = struct('rate', [], 'distribution', [], 'threshold', [], 'beamforming', []);

instanceSet = 1 : nInstances;
for iInstance = 1 : nInstances
	try
		Result(:, iInstance) = load(strcat(directory, 'instance_', num2str(iInstance)), 'Result').Result;
	catch
		Result(:, instanceSet == iInstance) = [];
		instanceSet(instanceSet == iInstance) = [];
	end
end
