function mainDOGS0(nit)
% Clear variables; be more agressive if the user has specified a preference to do so
% clear cleanupTask; clear % important to clear the cleanupTask variable before everything else
% clear all; close all; clc;

%% Load the "forecast" config file
%  This config file contains details such as
%	* Number of Prameters
%	* Optimization type (method, Number of Constraints)
%	* optVar: y0 or K
%	* bounds for x: bnd1, bnd2
%	* output directory


% User message
fprintf( 1 , 'Delta DOGS Z OMEGA OPTIMIZATION PACKAGE\n\n' );
fprintf( 1 , 'Configuring optimization:              ' ); timer_ = tic;

% functions can use a persistent variable to avoid repeatedly needing to re-read their conf file, but this means that we'll need to call 'clear functions' to force them to re-read it if we edit a conf file.  Which shouldn't be too bad?
[ conf, steps ] = DOGSConf( '../Zomega/DOGS.conf' );
global lattice n bnd1 bnd2 Search tri Ain bin
%%
n = str2num(conf.NumParam); ms = str2num(conf.NumConstraints); m=2*n; % initial starting points
bnd1 = str2num(conf.bnd1)'; bnd2 = str2num(conf.bnd2)';
Search.method = str2num(conf.method); Search.constant = str2num(conf.optVar);
iter_max = str2num(conf.MAX_ITER); MESH_SIZE = str2num(conf.MESH_SIZE);
% interpolation strategy:
inter_method = str2num(conf.interpolation_strategy);
% mesh grid points
lattice = conf.lattice;
MESHSIZE = str2num(conf.MESH_SIZE);

% load opt_prev
% nit = 0
if nit == 0
    %% calculate the initialization
    
    %find initial points vetices
    delta0=0.15;% how close the intial points needs to be from the body center
    [xU,xE,Ain, bin, length_scale,deltaU] = DOGS_initial_pts(bnd1,bnd2,n, delta0);
    %find closest mesh points to the initial pts
    xE = Unconstraint_quantizer(xE,MESHSIZE./length_scale);
    ninit = size(xE,2);
    neval=ninit;
    save pts_to_eval xE -ASCII
    save surr_x_eval_set xE -ASCII
%     save surr_J_eval_set yE -ASCII
    %     save surr_x_support_set xU -ASCII
    fid = fopen('store_neval','w');
    fprintf(fid,'%2i\n',neval);
    fclose(fid);
    stop_opt = 0;
    nit=1;
    save opt_prev.mat Search nit xE xU stop_opt MESHSIZE bnd1 bnd2 lattice
    fid = fopen('stop_file','w');
    fprintf(fid,'%2i\n',stop_opt);
    fclose(fid);
    return
else
    %read in new data points
    %has current best improved?
    clear yE xE xU
    ym = load('surr_J_new'); %store cost function result

    load('opt_prev')
    if nit==0,
        yE=[yE,ym];
          save surr_J_eval_set yE -ASCII
    end
    yE = load('surr_J_eval_set');
    y0=Search.constant;


    
    inter_par = interpolateparametarization(xE,yE,1);
    
    

    %%%
       
   % Calculate the metric of the unevaluated points.
    yup=[]; yrp=[];
     for ii=1:size(xU,2)
               yup(ii)=(interpolate_val(xU(:,ii),inter_par)-y0)/mindis(xU(:,ii),xE);
               yrp(ii)=interpolate_val(xU(:,ii),inter_par);
      end
          %  keyboard
   % check the minimum value of xU
            if size(xU,2)>0
                 [tup,ind_up]=min(yup);
            else
                 tup=inf;
            end
    
%% minimizing the search function
     xm=tringulation_search_bound(inter_par,[xE,xU],yE);
     
    %% 
    ninit = size(xE,2);
    neval=ninit;
    save pts_to_eval xm -ASCII
    save surr_x_eval_set xE -ASCII
    %     save surr_x_support_set xU -ASCII
    fid = fopen('store_neval','w');
    fprintf(fid,'%2i\n',neval);
    fclose(fid);
    stop_opt = 0;
    save opt_prev.mat Search nit xE xU stop_opt MESHSIZE bnd1 bnd2 lattice
    fid = fopen('stop_file','w');
    fprintf(fid,'%2i\n',stop_opt);
    fclose(fid);
    return
%%

end
    

end