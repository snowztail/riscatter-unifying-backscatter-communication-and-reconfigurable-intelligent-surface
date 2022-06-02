if ~exist(directory, 'dir')
	mkdir(directory)
end
if exist('iInstance', 'var')
	save(strcat(directory, 'instance_', num2str(iInstance)));
end
