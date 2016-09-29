%% Delta DOGS Omega
%
%  Title: Delta DOGS optimization configuration
%
%  Author: Shahrouz Alimohammadi
%
%  Description:
%    Gathers configuration code in one location to remove clutter from the
%    main optimization script
%
%
function [ conf , method] = DOGSConf( cfName )

% conf = readConf(DOGSGetConfPath( char(cfName) ),0);
conf = readConf(char(cfName) ,0);

if( ~isfield(conf,'NumParam') || ~isfield(conf,'method') || ~(isfield(conf,'optVar') || ~isfield(conf,'constraintType')) || ~isfield(conf,'NumConstraints') || ~isfield(conf,'outputDir') || ~isfield(conf, 'bnd1') || ~isfield(conf, 'bnd2'))
    % probably should break this out by field so we can be more specific
    error 'one of the required config parameters was not specified';
end
if( ~isfield(conf,'MESHSIZE') )
    conf.MESHSIZE = 8; % or some other suitable default value
    warning('DOGS:missingConfig_MESHSIZE','Initial mesh size threshold was not specified, so we are using the default (%i).  You should probably override this in DOGS.conf',conf.MESHSIZE);
else
    conf.MESHSIZE = str2num(conf.MESHSIZE); %#ok<ST2NM>
end
if( ~isfield(conf,'MAX_ITER') )
    conf.MAX_ITER = 500;
    warning('DOGS:missingConfig_MAX_ITER','No Maximum iteration was specified, so we will be doing ''%s''.  If this was not your intent, please specify a step in DOGS.conf',conf.MAX_ITER);
else
    conf.MESHSIZE = str2num(conf.MAX_ITER); %#ok<ST2NM>
end
if( ~isfield( conf , 'interpolaion_strategy' ) )
    conf.interpolaion_strategy = 'polyharmonic-spline';
    warning('DOGS:missingConfig_interpolaion_strategy','No interpolation strategy was specified, so we will be doing ''%s''.  If this was not your intent, please specify a step in DOGS.conf',conf.interpolaion_strategy);
end

% set a default value for the delta_tol param if it's missing. The criteria
% for stopping the algorithm
if( ~isfield(conf, 'delta_tol') )
    conf.delta_tol= 0.01;
    warning('DOGS:missingConfig_delta_tol','No minimum delta x was specified, so we will be doing ''%s''.  If this was not your intent, please specify a step in DOGS.conf',conf.delta_tol);
else
    conf.MESHSIZE = str2num(conf.delta_tol); %#ok<ST2NM>
end

method = strcat(conf.method, ' --- ', conf.constraintType);
end
