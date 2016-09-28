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
clc
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=2; m=2*n; 
% iter_max=10;
iter_max=30;
% iter_max=11;
RES = 1e-2;
% RES = 1e-3;
x_star=[0.7;0.1];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%%
Method = 2;
% Var = 0.1;
Var =0.08; % maaxiter 40
% Var =1;
%%
% Var = 0.2;
% Method =1;
%%
% Var = 0.4922;
%  Var = 2;
%  Var =10;
% Var = 1;
 % Var = 0.5;
% Var=0.3651;
%%
% Kf =0.75;
% Kf =8;
Kf=1;
Ex=2;
%%
if Ex==1
% end
%
% if nargin <3
% specify the objective function
% % % % % % % % % % % % % % % fun =@(x)x(2,:);
fun =@(x)(1-x(2,:)).*Kf;
%fun=@(x) (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2;
%%
% specify the contraints functions
b =3/4;
 con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
 con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%%
%con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%con{1}=@(x) (x(2,:)-con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%%
elseif Ex==2
KC2 = 1;
% KC2 = 1/3;
fun =@(x)x(2,:).*Kf;
con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
%%
end
    lob=zeros(n,1); upb=ones(n,1);
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

xi(1,1)=xi(1,1)+0.00001; % bc it couldn't calculate deluany triangulation
        for k=1:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
%   keyboard
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
           [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%
% %%%%%%%%%%%%%% FIGURE
%   figure(1); clf;
%     subplot(2,1,1)
%     plot(1:length(yi(1:end)),yi,'-')
%     title('objective function')
%     grid on
%     subplot(2,1,2)
%     plot(1:length(yi(1:end)),C{1},'-')
%     title('costrained function')
%     grid on
%     %
%    figure(2); clf;
% %    subplot(2,1,1)
%    tri=delaunayn(xi.');
%     triplot(tri,xi(1,:),xi(2,:))
%     hold on
    tt = 0:0.01:1;
%   %  keyboard
%    % plot(tt, fun_weierstrass(tt), 'r', 'linewidth', 3)  
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
% %    plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
% %     end
%     plot(xi(1,:),xi(2,:),'s') 
%     if Method ==1
%     title(strcat('Position of the points',' for K = ', num2str(Search.constant)))
%     elseif Method ==2
%    title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
%     end
%     xlim([0 1])
%     ylim([0 1])
% %      xlim([-2 2])
% %     ylim([-2 2])
%     drawnow   
% % Ff(k) = getframe;
% %     keyboard
%     
%    plot_figures_3(xi,yi,C, tt, Ex)
    %%
    
    
%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
if impro<RES|| k == iter_max+1

    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(33));
    % keyboard
%      save(strcat(fName,'_workspace','.mat')) 
%   csvwrite(strcat(fName,'.csv'), iter_1);
   break
end
%tt = 0:0.01:1;
%plot_figures_3(xi,yi,C, tt)
%%% interpolation plot at 0
  %xx=0:0.05:1;
%   for ii=1:length(xx)
%       yg1(ii)=interpolate_val([xx(ii);0.1],inter_par_g{1});
%       yg2(ii)=interpolate_val([xx(ii);0.1],inter_par_g{2});
%       yc1(ii)=con{2}([xx(ii);0.1]);
%   end
% figure(3)
% plot(xx,yg2,'--',xx,yc1,'-')

% figure(4)
% plot(1:k,y_star,'-',1:k,-1./ym,'--')
% drawnow

%x_prev = xm;
%y_prev = ym;
% keyboard

%%
clear U_f p_f G_c C_l Err TT SS G1 G2
% xv=0:0.01:1;
xv=0:0.02:1;
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
% % % %        CV(ii,jj) = max([G1, G2])+1./ym(k)* Err(ii,jj);
CV(ii,jj) = max([G1, G2]) - ym(k)* Err(ii,jj);
       if CC(ii,jj) >0
           CC(ii,jj) = nan;
       end
         if CV(ii,jj) >0
           CV(ii,jj) = nan;
       end
     
    end
end
%%
% figure(10); clf;
% % hold on
% % subplot(2,1,2)
% contourf((CC).')
% caxis([-0.5,0])
% colorbar
% title('g(x) - K e(x)')
% % Gf(k) = getframe;
% % pause
% 
% 
% figure(11)
% contourf((P_f-Var*Err).')
% caxis([-0.5,0])
% colorbar
% title('p(x)-Ke')
% 
% % pause
% % figure(3)
% % contourf(C_l.',-0.5:0.5:0)
% % % surf(C_l)
% % title('c(x)')
% % pause
% % figure(12)
% % % contourf(G_c.',-0.5:0.5:0)
% % contourf(C_l.')
% % title('c_l(x)')
% % pause
%%
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
% error plot in the domain 
% surf(Err.')
% contourf(Err.')
% title('err(x)')
% figure(14);clf;
	% Create axes
% 		ah(1) = axes('Position', [0.015 0.5 0.3 0.47]);
        		ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
contourf(Err.')
colorbar
title('e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
        ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
if Search.method ==1
%         contourf((CC).')
        %%
         hold on
%         contourf(xv,xv,(CC).', 0.1:0.1)
          contourf(xv,xv,(CC).')
  if Ex ==1
       plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
    elseif Ex==2
       plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
  end
     plot(xi(1,:),xi(2,:),'s') 
    hold off
    %%
caxis([-0.5,0])
colorbar
title('g(x) - K e(x)')
axis off

elseif Search.method ==2
    hold on
 contourf(xv,xv,(CV).')
caxis([-0.5,0])
colorbar
title('g(x) - K e(x)')
axis off   
 if Ex ==1
       plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
    elseif Ex==2
       plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
  end
     plot(xi(1,:),xi(2,:),'s') 
    hold off
%     
%    contourf(SS.')
% title('s(x) = - e(x) / T(x)')
%   colorbar   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
if Search.method ==1
contourf((P_f-Var*Err).')
caxis([0,1])
colorbar
title('p(x) - K e(x)')
axis off
elseif Search.method ==2
    %%    
% % % % contourf((P_f+1/ym(k)*Err).')
contourf((P_f - ym(k)*Err).')
% caxis([0,1])
colorbar
title('p(x) - K e(x)')
axis off
    %%
%  contourf(TT.')
% title('T(x) = max ( p(x) - y_0 , max(g_l(x)) )')
% caxis([0,1])
% colorbar
% axis off    
end
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
% % %         title(strcat('Position of the points',' for K_a = ', num2str(-1./ym(k))))
            title(strcat('Position of the points',' for K_a = ', num2str(ym(k))))

%    title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
		set(ah, 'FontSize', op.fontsize);
		drawnow
        
%        F(k) = getframe(h); 
%        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
%        imagerDir = './images/';
%        videoDir = '.'
% 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% 		encodeFramesToVideo(videoDir, 'png', mp4filename);


%%
disp( [' finished iteration ' ,num2str(k), ' ... '] )
% pause
        end
        
       
        