addpath(genpath(pwd));

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultaxesfontname', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');
set(groot, 'defaultlinelinewidth', 2);
set(groot, 'defaultstemlinewidth', 2);
set(groot, 'defaultscatterlinewidth', 2);

% PBS_ARRAY_INDEX is environment variable: 0 means local and positive integer means HPC instance index
iInstance = str2double(getenv('PBS_ARRAY_INDEX'));
disp(iInstance);
% rng(iInstance);
rng shuffle;
