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
clear all; close all; clc;
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=2; m=2*n; 
%iter_max=10;
 iter_max=20;
RES = 1e-2;
x_star=[0.7;0.1];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%%
Method = 2;
Var = 0.2;
% % Method =1;
% %  Var = 0.3;
% Var = 1;
% Var = 1;
% Kf =0.75;
% Kf = 100;
% Kf = 2;
Kf=1;

% Kf =2;
Ex=1;
%%
if Ex==1
% end
x_star = [0.5,0.8];
%
% if nargin <3
% specify the objective function
% % % % % % % % % % % % % % % fun =@(x)x(2,:);
fun =@(x)(1-x(2,:)).*Kf;
%fun=@(x) (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2;
%%
% specify the contraints functions
% b =3/4;
% b =2/4;
% b =0.5;
b =0.8;
%  con{1}=@(x)x(2,:).^2+sin(x(1,:).*pi/2*2).^2;
%  con{2}=-con{1};
 con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
 con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%%
elseif Ex==2
% KC2 = 2;
% % KC2 = 1/3;
fun =@(x)x(2,:);
con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
end
%%   control constraints
    lob=zeros(n,1);
    upb=ones(n,1);
%%
Search.method = Method;
Search.constant = Var;
bnd1 = lob;
bnd2 = upb;
% Search.method=1;
% Search.constant=1;
% interpolaion strategy
inter_method=1;
% upper and lower bouds

% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
  xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
for ii=1:m
    acon{ii}=[];
end
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
delta_tol=0.2;

xi(1,1)=xi(1,1)+0.00001;
        for k=2:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
%             [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%
% %%%%%%%%%%%%%% FIGURE
    tt = 0:0.01:1;
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

   break
end


%%
clear U_f p_f G_c C_l Err TT SS G1 G2
xv=0:0.01:1; 
for ii=1:length(xv)
    for jj=1:length(xv)
%         U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
         P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
         C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
%         G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
        G1=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
        G2=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{2});
       [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
       CC(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
%          CC_2(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
% %        if CC(ii,jj) >0
% % %            CC(ii,jj) = nan;
% % CC(ii,jj) = 0;
% %        else 
% %            CC(ii,jj) = 1;
% %        end
       
       
    end
end
%%

h=figure; clf;
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
%%
% ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
% figure(10); clf;
axis square
box on

       hold on
%  contourf(xv,xv,(CC).', [1, 0.00000000001], 'linestyle', 'none')
% contourf(xv,xv,A);
%  contourf(xv,xv,CC.', [-10:10:0]);
contourf(xv,xv,-CC.', 0:10:10,  'linestyle', 'none' );
% % % contourf(xv,xv,(CC).', 0:1:1, 'linestyle', 'none')
  grid on
  grid minor
   plot(x_star(1),x_star(2),'kp','MarkerFaceColor','b', 'MarkerSize', 18)  
brighten(0.5)
% colormap(map)
colormap('bone')
% colorbar

% caxis([-.010,0])
% contourf(xv,xv,(CC).')
% contourf(xv,xv,(CC).','linestyle', 'none',0.1:0.1)
% set(h,'linestyle','none');
% caxis([-0.0001,0.0001])
% colormap(map)
  %%
%   plot(xi(1,:),xi(2,:),'ro', 'MarkerSize', 13) 
 if Ex ==1
       plot(tt, b*sin(pi*tt), 'k', 'linewidth', 3)
    elseif Ex==2
       plot(tt, rastriginn(tt), 'k', 'linewidth', 3)  
  end
  %% plot the evaluated points
  plot(xi(1,1:4),xi(2,1:4),'ks',  'MarkerFaceColor','k','MarkerSize', 15) 
  plot(xi(1,5:end),xi(2,5:end),'ks','MarkerFaceColor','k', 'MarkerSize', 15)
%   plot(xi(1,4:end-1),xi(2,4:end-1),'ro', 'MarkerSize', 13) 
    plot(xi(1,end),xi(2,end),'wx', 'MarkerSize', 15) 
  %%
 
 %%
   
% caxis([-0.5,0])
% colorbar
% title(['iter = ', num2str(k)], 'FontSize', 20)
% title('max (g_l(x) - K e(x))')
% axis off

set(gca, 'ytick', [])
set(gca, 'xtick', [])
%  set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 18);
% saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
% saveas(h,strcat('./images/', num2str(k,'%02d')), 'fig')
% saveas(h,strcat('./images/', num2str(k,'%02d')), 'epsc2')
% cmfile = strcat('./images/', num2str(k,'%02d'));
% figure_to_publish(cmfile)





%%
% % 
% % h2=figure(14); clf;
% % set(h2,'units','normalized','outerposition',[0 0 0.9 0.9])
% % %default values
% % op.fontsize = 12;
% % op.quiet = 0;
% % % op.dirpath = conf.outputDir;
% % op.dirpath = './';
% % op.width = 1500;
% % op.height = 900;
% % op.visible = 1;
% % op.title = 1;
% % op.cache = 0;
% % op.output = 0;
% % op.embed = 0;
% % op.save = 0;
% % op.label = 1; %Settign to zero removes labels on figures
% % %%
% % % ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
% %        hold on
% % contourf(xv,xv,(P_f-Var*Err).')
% % colorbar 
% % caxis([-0.1,1.])
% % 
% % % caxis([-0.5,0])
% % % colorbar
% % % title(['iter = ', num2str(k)], 'FontSize', 20)
% % % title('max (g_l(x) - K e(x))')
% % axis off
% % % saveas(h,strcat('./images_p/', num2str(k,'%02d')), 'png')
% % % saveas(h,strcat('./images_p/', num2str(k,'%02d')), 'fig')
% % % 
% % % 

drawnow
% %%
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
% 
% ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
% contourf(Err.')
% colorbar
% title('e(x)')
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
%         ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
%         hold on
%   contourf((CC).')
%   %%
%   plot(xi(1,:),xi(2,:),'rs', 'MarkerSize', 13) 
%   if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
%     elseif Ex==2
%        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
%   end
%  %%
% caxis([-0.5,0])
% colorbar
% title('g(x) - K e(x)')
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% % ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
% ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
% contourf((P_f-Var*Err).')
% caxis([0,1])
% colorbar
% title('p(x) - K e(x)')
% axis off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
%         	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
%     plot(1:length(yi(1:end)),yi,'-')
%     title('objective function')
%     grid on
%   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
%         	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
%         hold on
%     triplot(tri,xi(1,:),xi(2,:)) 
%     if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
%     elseif Ex==2
%        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
%     end
%     plot(xi(1,:),xi(2,:),'s') 
%     hold off
%       if Method ==1
%     title(strcat('Position of the points',' for K = ', num2str(Search.constant)))
%     elseif Method ==2
%    title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
%     end
%     xlim([0 1])
%     ylim([0 1])
% 		set(ah, 'FontSize', op.fontsize);
% 		
%         
% %        F(k) = getframe(h); 
%        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
% %        imagerDir = './images/';
% %        videoDir = '.'
% % 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% % 		encodeFramesToVideo(videoDir, 'png', mp4filename);
%%
% % % h=figure(13); clf;
% % % set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
% % % %default values
% % % op.fontsize = 12;
% % % op.quiet = 0;
% % % % op.dirpath = conf.outputDir;
% % % op.dirpath = './';
% % % op.width = 1500;
% % % op.height = 900;
% % % op.visible = 1;
% % % op.title = 1;
% % % op.cache = 0;
% % % op.output = 0;
% % % op.embed = 0;
% % % op.save = 0;
% % % op.label = 1; %Settign to zero removes labels on figures
% % % 
% % %         		ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
% % % contourf(Err.')
% % % colorbar
% % % title('e(x)')
% % % axis off
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
% % %         ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
% % %   contourf((CC).')
% % % caxis([-0.5,0])
% % % colorbar
% % % title('g(x) - K e(x)')
% % % axis off
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% % % % ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
% % % ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
% % % contourf((P_f-Var*Err).')
% % % caxis([0,1])
% % % colorbar
% % % title('p(x) - K e(x)')
% % % axis off
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
% % %         	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
% % %     plot(1:length(yi(1:end)),yi,'-')
% % %     title('objective function')
% % %     grid on
% % %   
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
% % %         	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
% % %         hold on
% % %     triplot(tri,xi(1,:),xi(2,:)) 
% % %     if Ex ==1
% % %        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
% % %     elseif Ex==2
% % %        plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
% % %     end
% % %     plot(xi(1,:),xi(2,:),'s') 
% % %     hold off
% % %       if Method ==1
% % %     title(strcat('Position of the points',' for K = ', num2str(Search.constant)))
% % %     elseif Method ==2
% % %    title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
% % %     end
% % %     xlim([0 1])
% % %     ylim([0 1])
% % % 		set(ah, 'FontSize', op.fontsize);
% % % 		
% % %         
% % % %        F(k) = getframe(h); 
% % %        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
% % % %        imagerDir = './images/';
% % % %        videoDir = '.'
% % % % 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% % % % 		encodeFramesToVideo(videoDir, 'png', mp4filename);
% pause(1)
% % keyboard
% % figure(10),clf;
% % C4= CC>0;
% % % C4(CC<=0)=0.1;
% % C3 = ~C4;
% % imagesc(flipud(C3.'))
% % imshow(flipud(C3.'))

  grid on
  grid minor
  view(0,-90)
% Path = 'C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\ex0\';
Path = 'C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\paper_nonlinear\ex0\';
% save(strcat('ex0_simulation_results_iter_', num2str(k), '.mat'))
% saveas(gcf,strcat(Path, num2str(k-1, '%02d')), 'fig')
% saveas(gcf,strcat(Path, num2str(k-1, '%02d')), 'epsc2')
% saveas(gcf,strcat(Path, num2str(k-1, '%02d')), 'png')
%%
        end
        
       
        