% Clear variables; be more agressive if the user has specified a preference to do so
clear cleanupTask; clear % important to clear the cleanupTask variable before everything else
clear all; close all; clc;

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
[ conf, steps ] = DOGSConf( 'DOGS.conf' );
global lattice n bnd1 bnd2 Search tri
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

load opt_prev
stop_opt=0;
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
    yE = load('surr_J_new'); %store cost function result
    load('surr_C_new'); %store cost function result
    
    load('opt_prev')
    y0=Search.constant;
    
    
    for k=1:iter_max
        
        
        % Modify the interpolations
        inter_par_p= interpolateparametarization(xE,yE,inter_method);
        for jj=1:ms
            inter_par_g{jj}= interpolateparametarization(xE,C{jj},inter_method);
        end
        
        
        % Calculate Discrete search function
        yup=zeros(1,size(xU,2));
        
        
        
        for ii=1:size(xU,2)
            %yup(ii)=estimate_max_cons_val(xU(:,ii),inter_par_p,inter_par_g,y0,ms)/mindis(xU(:,ii),xE);
            yup(ii)=inf;
        end
        
        
        
        % Perform the search
        while 1
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,[xE xU]);
            %keyboard
            xm = (xm-bnd1)./(bnd2-bnd1).*MESH_SIZE; xm = Unconstraint_quantizer(xm,1);
            xm =xm./ MESH_SIZE.*(bnd2-bnd1)+bnd1;
            %     xm=round((xm-bnd1)./(bnd2-bnd1).*MESH_SIZE)./MESH_SIZE.*(bnd2-bnd1)+bnd1;
            [xm,xE,xU,newadd,success]=points_neighbers_find(xm,xE,xU);
            
            
            if success==1
                break
            else
                yup=[yup estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)];
            end
            
            
        end
        
        
        if (estimate_max_cons_val(xm,inter_par_p,inter_par_g,y0,ms)/mindis(xm,xE)>min(yup) || mindis(xm,xU)<1e-6)
            [t,ind]=min(yup);
            xm=xU(:,ind); xU(:,ind)=[];
        end
        
        if mindis(xm,[xE xU])<1e-6
            save opt_prev.mat
            fid = fopen('stop_file','w'); fprintf(fid,'%2i\n',stop_opt); fclose(fid);
            return;
        end
        
        save pts_to_eval xm -ASCII
        
        %% TODO COMPLETE THE REST OF THE CODE......
        %[con , FUN, CS , CON, Xm] = IMEXRK_Solver3(xm',DX);
        %% Stoping Criteria
        %impro_p=mindis(xmp,xE);
        %impro_q=mindis(xmq,xE);
        %
        [ViolCon ] = const_violation( con, ConBound );
        %
        conV = cell2mat(ViolCon)>0;
        CX = [con sqrt(FUN.val)];
        %     keyboard
        ConStop = max(0, max(cell2mat(ViolCon)));
        %  if impro<RES|| k == iter_max +1 || ConStop ==0
        xE=[xE xm];
        %if isx ~=1
        %   xU=[xU xm];
        %end
        yE = [yE FUN.val];
        for jj=1:ms
            C{jj}=[C{jj} real(con{jj})];
        end
        % figure(1)
        % plot(yi)
        % drawnow
        %% ploting constraints
        % plot_preview(xi,yi,C )
        %     prettyplot_IMEXRK( xi, yi, C )
        [  CV_S ] = const_violation( C, ConBound );
        prettyplot_IMEXRK( xE, yE, CV_S)
        drawnow
        disp([num2str(k/iter_max*100), ' % Completed'])
    end
    
end
% User message
fprintf( 1 , '%5.2f seconds.\n' , toc(timer_) );
