
%global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
global m
conf_path = './DOGS.conf';
conf = readConf('./DOGS.conf',0);

%
opt_init()
%
RES = str2num(conf.delta_tol);
iter_max = str2num(conf.MAX_ITER);
m = str2num(conf.NumParam)*2;
%
nit = 1;
if nit ==0
 %find initial points using LHS
    [xi, acon] = initial_pts();
    xi=xi';
    %find closest mesh points to the LHS initial pts

    neval=m;
    save pts_to_eval xi -ASCII
    fid = fopen('store_neval','w');
    fprintf(fid,'%2i\n',neval);
    fclose(fid);
   prev_it =1;
    cost_improve=1;
    curr_bestJ = 1000000;
    xm = [];
    stop_opt = 0;
    nit=1;
    save opt_prev nit  cost_improve curr_bestJ xm stop_opt 
    fid = fopen('stop_file','w');
    fprintf(fid,'%2i\n',stop_opt);
    fclose(fid);
%     return
keyboard
    %
    
    %
else
  %%  
 %read in new data points
 %has current best improved?
 ym = load('surr_yi_new') %store cost function result
 xm = load('surr_xi_new') %Design parameter
    load surr_pts
    improve = 0;
    Jsurr = [];
    Xsurr = [];
    neval = size(ym,1);
    for i = 1:length(ym)
        if yi(i) ~= 1000000
            Jsurr = [Jsurr; ym(i)];
            Xsurr = [Xsurr; xm(i,:)];
        end
        if yi(i) < curr_bestJ
            curr_bestJ = yi(i)
            curr_bestX = xi(i,:)
            disp('update current best point')
            curr = [curr_bestJ curr_bestX nit neval];
            if nit == 1
                save best_Jhist curr -ASCII
            else
                save best_Jhist curr -ASCII -append
            end
            cost_improve = 1;
            improve = 1;
        end
    end
    keyboard
      if improve == 0 & nit > 1
        cost_improve = 0;
    elseif nit == 1
        cost_improve = 1;
    end
    % save optimization history 
    if nit == 1
        save it_hist prev_it -ASCII
    else
       size_prepts=size(xi);
        it_state=[prev_it,size_prepts(1),cost_improve]; %record the type, numbers and result of last evaluations  
        save it_hist it_state -ASCII -append
    end
      yi = [yi; Jsurr];
    xi = [xi; Xsurr];
    save surr_pts yi xi
    nit = nit+1;
end
    
[yi,xi] = dsmerge(xi, yi);
x_prev = xi(:,1);
y_prev = yi(:,1);




 for k=2:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
%             [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%
%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
if impro<RES|| k == iter_max +1

    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(33));
    % keyboard
%      save(strcat(fName,'_workspace','.mat')) 
 fid = fopen('stop_file','w');
    fprintf(fid,'%2i\n',1);
 fclose(fid);
   return
end
    
    
 end

