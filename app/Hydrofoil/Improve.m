function [y,sigma] = Improve(ind,t)

% Improve the measurement at Data point with index ind

% go to the Folder

% run the NFA: Modify the restart file, input file, time series

% Read the time series file
[t,zf]=Read_log_file('Zforces.log');

% statical analysis calculation
Initial_mod(1)=1; Initial_mod(2)=t;
[mu,sigma,Initial_mod] = alpha_statitical_analyze(zf,Initial_mod);
% Improve the time scale.

    
end





