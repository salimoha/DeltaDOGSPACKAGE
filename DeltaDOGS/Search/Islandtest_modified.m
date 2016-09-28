% 2d Islands test problem
% this is for adoptive K algorithm
clear all
close all
clc
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=1; m=2*n; b =1/2;
% iter_max=26;
iter_max=36;
% iter_max=60;
RES = 1e-2;
x_star=[0.7;0.1];
x0=[0;0];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%%
Method = 2;
Var = 0.2;
% Var = 0.6;
% % % % Method = 1; %TODO
% % % % Var =1;
%%
% Var = 0.3;
% KC2 = 10;
KC2 = 1;
KCf = 1;
% Method = 1;
% % Var = 0.2;
% Var = 1;
% end
%
% if nargin <3
% specify the objective function
%fun =@(x)x(2,:);
fun=@(x) ( (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2 )*KCf;
% specify the contraints functions
% con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
% con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%con{1}=@(x) (x(2,:)-con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));

%con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
%con{2}= @(x) -con{1}(x);
% con{1}=@(x) rastriginn2(x)-0.5;
% con{1}=@(x)fun_constraint;
con{1}=@(x) (rastriginn2(x)-0.5 )*KC2; %the same
%fun=@fit8;
%fun=@plygen;
%fun=@rastriginn;
% end

% if nargin <5
    lob=zeros(n,1); upb=ones(n,1);
%     lob=zeros(n,1); upb= [ 1.7;2];
% end
%
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
acon =cell(1,m);
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
        for k=1:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%%%%%%%%%%%%% FIGURE
% % % 
% % % % % % % %   figure(1); clf;
% % % % % % % %     subplot(2,1,1)
% % % % % % % %     plot(1:length(yi(1:end)),yi,'-')
% % % % % % % %     title('objective function')
% % % % % % % %     grid on
% % % % % % % %     subplot(2,1,2)
% % % % % % % %     plot(1:length(yi(1:end)),C{1},'-')
% % % % % % % %     title('costrained function')
% % % % % % % %     grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   h=figure(2); 
%    clf;
   tri=delaunayn(xi.');
%     triplot(tri,xi(1,:),xi(2,:))
    hold on
   tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.5:0.5)
colormap gray
% colormap bone
brighten(0.5)
%     end
strmax = ['   ', num2str(k)];
% text(xmax,ymax,strmax,'HorizontalAlignment','right');
text(xi(1,end),xi(2,end),strmax ,'HorizontalAlignment','right');
    plot(xi(1,end),xi(2,end),'rs') 
    if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant) ))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
%      xlim([-2 2])
%     ylim([-2 2])
   triplot(tri,xi(1,:),xi(2,:))
drawnow   
% % % 
% % % %     impro=norm(x_prev-xm);
% % %       impro=mindis(xm,xi(:,1:end-1));
% % % %   saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
tt = 0:0.01:1;
% % % plot_figures_3(xi,yi,C, tt)
%% interpolation plot at 0
%   xx=0:0.05:1;
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



%% movie maker
xv=0:0.1:1; 

for ii=1:length(xv)
    for jj=1:length(xv)
        U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
        P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
        C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
        G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
       [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
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
% %          contourf((P_f-Var).')
% % title('p(x) - y_0')
contourf(G_c.')
title('g(x)')
% caxis([-0.5,0])
colorbar
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);

contourf(TT.')
title('T(x) = max ( p(x) - y_0 , max(g_l(x)) )')
% caxis([0,1])
colorbar
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
        	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
contourf(SS.')
title('s(x) = - e(x) / T(x)')
  colorbar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
        	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
       hold on
            tri=delaunayn(xi.');
    triplot(tri,xi(1,:),xi(2,:))

   tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.5:0.5)
colormap gray
% colormap bone
brighten(0.5)
%     end
    plot(xi(1,:),xi(2,:),'rs') 
    xlim([0 1])
    ylim([0 1])
   triplot(tri,xi(1,:),xi(2,:))
drawnow   

%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
		set(ah, 'FontSize', op.fontsize);
		        
%        F(k) = getframe(h); 
%        saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
%        imagerDir = './images/';
%        videoDir = '.'
% 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% 		encodeFramesToVideo(videoDir, 'png', mp4filename);












        end
        
        
% % %         
% % % rootpath = './figures_paper_nonlinear_island/';
% % % if ~exist(rootpath) == 1 
% % % mkdir(rootpath);
% % % end
% % % EX=2;
% % % Var = Search.constant;
% % % Method = Search.method;
% % % 
% % %         
% % % fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
% % % if ~exist(fpath) == 1 
% % % mkdir(fpath);
% % % end
% % % if Var <0
% % %         % file figure to be saved
% % %     fName = strcat( fpath ,'/Method_',num2str(Method), '_y0_neg_', num2str(abs(Var)*100),'KC2_' ,num2str(KC2), '_ex_' ,num2str(EX));
% % %     % file figure to be saved
% % % else
% % %     fName = strcat( fpath ,'/Method_',num2str(Method), '_y0_', num2str(Var*100),'KC2_' ,num2str(KC2), '_ex_' ,num2str(EX));
% % %     
% % % end
% %     consraints_delta_dogs_3(Var, Method, fun, con);
%     %%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% %     figure_to_publish(fName, strcat('Method ', num2str(Method), ' K = ', num2str(Var)  ))
%        figure_to_publish(fName)
% figure(1)
%     figure_to_publish(strcat(fName, '_convegance')) 
%         
        