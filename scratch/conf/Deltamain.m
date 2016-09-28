% Clear variables; be more agressive if the user has specified a preference to do so
clear cleanupTask; clear % important to clear the cleanupTask variable before everything else
clear all; close all; clc;
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
%% Load the "forecast" config file
%  This config file contains details such as
%	* start time
%	* stop time
%	* which step to run
%	* imager name
%	* deployment site name (or some other way to identify the entity that the forecast gets generated for)
%	* output directory

% User message
fprintf( 1 , 'Delta DOGS Z OMEGA OPTIMIZATION PACKAGE\n\n' );
fprintf( 1 , 'Configuring optimization:              ' ); timer_ = tic;
% functions can use a persistent variable to avoid repeatedly needing to re-read their conf file, but this means that we'll need to call 'clear functions' to force them to re-read it if we edit a conf file.  Which shouldn't be too bad?
[ conf, scheme ] = DOGSConf( 'DOGS.conf' );

% User message
fprintf( 1 , '%5.2f seconds.\n' , toc(timer_) );


outputDir = [conf.outputDir  '/' datestr(now,'yyyymmdd') '/'];
% Fix path for windows
if( ispc() ), outputDir = char(java.io.File(outputDir)); end

fprintf( 1 , '%5.2f seconds.\n' , toc(timer_) );
%
n = str2num(conf.NumParam); ms = str2num(conf.NumConstraints); m=2*n; % initial starting points
bnd1 = str2num(conf.amin)'; bnd2 = str2num(conf.amax)';
Search.method = str2num(conf.method); Search.constant = str2num(conf.optVar);
iter_max = str2num(conf.MAX_ITER); MESH_SIZE = str2num(conf.MESH_SIZE);
%%
[xU,xE,Ain, bin] = DOGS_initial_pts(bnd1,bnd2,n)

%% Initialization
% interpolaion strategy
inter_method=1;
% upper and lower bouds
 bnd1=zeros(n,1);
  bnd2=ones(n,1);
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
