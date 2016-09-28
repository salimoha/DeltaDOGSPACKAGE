% function [k] = consraints_delta_dogs(Var, Method, fun, con,lob,upb)
% The Delta_Dogs algorithm with nonlinear costraints.
% Simple bound costraints for the control parameter is imposed.
% bnd1<x<bnd2
% fun: The objective function
% con{i}= The costraint function
% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints
% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K
clear all
close all
% clc
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=2; m=2*n; 
%iter_max=10;
iter_max=75;
RES = 1e-3;
%  RES = 0.05;
x_star=[0.7;0.1];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%%
% % % % Method = 2;
% % % % Var = 0.1;
Method =1;
%  Var = 0.3;
% Var = 1;
% Var = 1;
% Kf =0.75;
% Kf = 100;
% Kf = 2;
Kf=1;
% Var = 0.1*Kf;
Ex=2;
for Method = 1
for Var =  0.3651
   
%%
if Method ==2
 Var = Var*Kf;
end
 % figure(gcf); clf;
clear xi yi inter_par_p inter_par_g tri C
if Ex==1
fun =@(x)(1-x(2,:)).*Kf;
b =0.8;
con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%%
elseif Ex==2
 KC2 = 1;
% % KC2 = 1/3;
fun =@(x)x(2,:)*Kf;
con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
end
%%   control constraints
lob=zeros(n,1);
upb=ones(n,1);
%% configuration for the method of constant K or adaptive K
Search.method = Method;
Search.constant = Var;
bnd1 = lob;
bnd2 = upb;
% Search.method=1;
% Search.constant=1;
% interpolaion strategy
inter_method=1;
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
acon = cell(1,m);
for ii=1:size(xi,2)
    ind=1:2*n;
    ind=ind(Ain*xi(:,ii)-bin>-0.01);
    for jj=1:length(ind)
    acon{ind(jj)}=[acon{ind(jj)} ii];
    end
yi(ii)=fun(xi(:,ii));
    for jj=1:ms
        C{jj}(ii)=con{jj}(xi(:,ii));
    end
end
x_prev = xi(:,1);
y_prev = yi(:,1);
delta_tol=0.1;
xi(1,1)=xi(1,1)+0.00001;
iter_max = 1;
        for k=1:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
% % % % % % % % % % %            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%
% %%%%%%%%%%%%%% FIGURE
    tt = 0:0.1:1;
%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
if impro<RES|| k == iter_max+1
    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures_paper';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(k,'%02d'));
    % keyboard
%     if ~exist(fpath) 
%      mkdir(fpath)
%     else
%         disp('this simulation is already made  ....\n')
%         break
%         
%     end
% % %  set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 18);
% saveas(h,fName, 'png')
% saveas(h,fName, 'fig')
% saveas(h,fName, 'epsc2')
%      save(strcat(fName,'_workspace','.mat')) 
   break
end
clear U_f p_f G_c C_l Err TT SS G1 G2
% % xv=0:0.01:1; 
% % for ii=1:length(xv)
% %     for jj=1:length(xv)
% %         U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
% %          P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
% %          C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
% %         G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
% %         G1=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
% %         G2=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{2});
% %        [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
% %        CC(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
% %          CC_2(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
% %        if CC(ii,jj) >0
% %            CC(ii,jj) = nan;
% % CC(ii,jj) = 0;
% %        else 
% %            CC(ii,jj) = 1;
% %        end
% %        
% %        
% %     end
% % end
% % 
if Search.method ==2
% % 
% %% movie maker
% xv=0:0.05:1; 
% for ii=1:length(xv)
%     for jj=1:length(xv)
%         U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
%         P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
%         C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
%         G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
%        [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
%     end
% end
% 
% %%
% h=figure(11); 
% % set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
% %default values
% op.fontsize = 12;
% op.quiet = 0;
% % op.dirpath = conf.outputDir;
% op.dirpath = './';
% op.width = 1000;
% op.height = 1000;
% op.visible = 1;
% op.title = 0;
% op.cache = 0;
% op.output = 0;
% op.embed = 0;
% op.save = 0;
% op.label = 1; %Settign to zero removes labels on figures
% hold on
% %   plot(xi(1,1:4),xi(2,1:4),'ks',  'MarkerFaceColor','g','MarkerSize', 16) 
% plot(xi(1,1:4),xi(2,1:4),'ks','MarkerSize', 16) 
%   %%
%   if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
% %         rectangle('Position',[0 0 1 1])      
%     elseif Ex==2
%        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
%   end 
% %% plot the evaluated points
% %   plot(xi(1,end),xi(2,end),'ks',  'MarkerFaceColor','g','MarkerSize', 16) 
% plot(xi(1,end),xi(2,end),'ks','MarkerSize', 16) 
% %%
% %  triplot(tri,xi(1,:),xi(2,:)) 
% strmax = ['   ', num2str(k)];
% % text(xmax,ymax,strmax,'HorizontalAlignment','right');
% text(xi(1,end),xi(2,end),strmax ,'HorizontalAlignment','right');
% box on
% % set(gca, 'ytick', [])
% % set(gca, 'xtick', [])
% % % cmfile = strcat('./images/', num2str(k,'%02d'));
% % % figure_to_publish(cmfile)
% drawnow
% h=figure(13); clf;
% set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
% %default values
% op.fontsize = 12;
% op.quiet = 0;
% % op.dirpath = conf.outputDir;
% op.dirpath = './';
% op.width = 1500;
% op.height = 900;
% op.visible = 1;
% op.title = 1;
% op.cache = 0;
% op.output = 0;
% op.embed = 0;
% op.save = 0;
% op.label = 1; %Settign to zero removes labels on figures
% % error plot in the domain 
% % surf(Err.')
% % contourf(Err.')
% % title('err(x)')
% % figure(14);clf;
% 	% Create axes
% % 		ah(1) = axes('Position', [0.015 0.5 0.3 0.47]);
%         		ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
% contourf(Err.')
% colorbar
% title('e(x)')
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
%         ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
% % %          contourf((P_f-Var).')
% % % title('p(x) - y_0')
% contourf(G_c.')
% title('g(x)')
% caxis([-0.5,0])
% colorbar
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% % ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
% ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
% 
% contourf(xv,xv,TT.')
% title('T(x) = max ( p(x) - y_0 , max(g_l(x)) )')
% % caxis([0,1])
% colorbar
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
%         	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
%             hold on
% contourf(xv,xv,SS.')
% title('s(x) = - e(x) / T(x)')
%   caxis([-1,2])
%       plot(xi(1,:),xi(2,:),'rs') 
%            plot(xi(1,end),xi(2,end),'kx') 
%     xlim([0 1])
%     ylim([0 1])
%     colorbar
%     hold off
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
%         	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
%        hold on
%             tri=delaunayn(xi.');
%     triplot(tri,xi(1,:),xi(2,:))
% 
%    tt=0:0.01:1;
% 
% 
%   %%
%   if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
% %         rectangle('Position',[0 0 1 1])      
%     elseif Ex==2
%        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
%   end 
%   
% % colormap bone
% % brighten(0.5)
% %     end
%     plot(xi(1,:),xi(2,:),'rs') 
%     xlim([0 1])
%     ylim([0 1])
%    triplot(tri,xi(1,:),xi(2,:))
%     plot(xi(1,end),xi(2,end),'kx') 
% drawnow   
% % pause
% % 
% % %     impro=norm(x_prev-xm);
% %       impro=mindis(xm,xi(:,1:end-1));
% % 		set(ah, 'FontSize', op.fontsize);
% % 		        
% % %        F(k) = getframe(h); 
% % %        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
% % %        imagerDir = './images/';
% % %        videoDir = '.'
% % % 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% % % 		encodeFramesToVideo(videoDir, 'png', mp4filename);
% % % 
% % % if k > 15
% % %   keyboard
% % % end
% % figure(4)
% % x1 = linspace(0,1,1000);
% % for ii=1:1000
% % y1(ii)=interpolate_val([ x1(ii);0.1],inter_par_g{1});
% % end
% %  plot(x1,y1,'linewidth',3)

        end
%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%
% CONSTANT K PLOTS
%%
if Search.method ==1
    
xv=0:0.02:1; 
for ii=1:length(xv)
    for jj=1:length(xv)
        U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
         P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
         C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
        G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
        G1=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
        G2=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{2});
       [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
       CC(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
         CC_2(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
       if CC(ii,jj) >0
           CC(ii,jj) = nan;
CC(ii,jj) = 0;
       else 
           CC(ii,jj) = 1;
       end
       
       
    end
end

    
    
h=figure(13); clf;
set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
%default values
op.fontsize = 12;
op.quiet = 0;
% op.dirpath = conf.outputDir;
op.dirpath = './';
op.width = 1500;
op.height = 900;
op.visible = 1;
op.title = 1;
op.cache = 0;
op.output = 0;
op.embed = 0;
op.save = 0;
op.label = 1; %Settign to zero removes labels on figures

ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
contourf(Err.')
colorbar
title('e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
        ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
        hold on
  contourf(xv,xv,(CC).')
  %%
  plot(xi(1,:),xi(2,:),'rs', 'MarkerSize', 13) 
  xm1 =[0.8852;0.1926];
  plot(xm1(1),xm1(2),'kx', 'MarkerSize', 13) 
  axis square
%   if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
%     elseif Ex==2
%        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
%   end
 %%
% caxis([-0.5,0])
colorbar
title('g(x) - K e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
contourf((P_f-Var*Err).')
caxis([0,1])
colorbar
title('p(x) - K e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
        	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
    plot(1:length(yi(1:end)),yi,'-')
    title('objective function')
    grid on
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
        	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
        hold on
    triplot(tri,xi(1,:),xi(2,:)) 
    if Ex ==1
       plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
    elseif Ex==2
       plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
    end
    plot(xi(1,:),xi(2,:),'s') 
    hold off
      if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant)))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
		set(ah, 'FontSize', op.fontsize);
		
        
%        F(k) = getframe(h); 
%        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
%        imagerDir = './images/';
%        videoDir = '.'
% 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% 		encodeFramesToVideo(videoDir, 'png', mp4filename);
    
end

        end
end
end

        