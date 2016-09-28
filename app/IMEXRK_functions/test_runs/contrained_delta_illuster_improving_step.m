% Representation of subinterpolation
clear all
close all
clc
warning('off','all')
f0 = 0.5455;
savekon=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 type='exploring'
type='replacing'
switch type
    case 'exploring'
        xi=[0.18 0.3 0.8]; % exploring step
    case 'replacing'
        xi=[0.18 0.5 0.8]; % replacing step
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c =-.2;
w=10;
fun=@(x) (1*(x-c).^2+sin(w*((x-c)+0.2)))/50;
% fun=@(x) (1*(x-c).^2+0.1*cos(w*((x-c)+0.2)))/20;
gfun=@(x) (1*(x-c).^2+1.*sin(w*((x-c)-0.1)))/20;
% xi=[0.001 1 0.2 0.3 0.5 0.6 0.12 0.07 .17 .15];
% xi=[0.2 0.7 0.5]; % replacing step
% xi=[0.2 0.3 0.5]; % replacing step
% xi=[0.1 0.7 0.3];
%  xi=[0.2 0.8 0.1];% exploring step

xiu=[0 1];
for ii=1:length(xi)
    yi(ii)=fun(xi(:,ii));
    gi(ii)=gfun(xi(:,ii));
end
%T=ones(1,3);
%K=2; L=1;
%[xi,ind]=sort(xi); yi=yi(ind);
xx=0:0.01:1;
for ii=1:length(xx)
    yy(ii)=fun(xx(ii));
    gg(ii)=gfun(xx(ii));
end
%%
figure(1); clf;
hold on
plot(xx,yy,'k--',xx,gg,'k-','linewidth',2)
% plot(xx,max(yy,gg)+1.,'k-.','linewidth',3)
plot(xx,max(yy,gg)+.0003,'-.','linewidth',4)
set(gca,'YTick',[])
set(gca,'XTick',[])
% axis square
% legend('f(x)','g(x)','max\{f,g\}')
% grid minor
xiT=[xi xiu]; xiT=sort(xiT);
% plot(xiu,max(fun(xiu),gfun(xiu)),'k*', 'Markersize',10, 'MarkerFacecolor', 'k')
% plot(xiu,min(fun(xiu),gfun(xiu)),'k*', 'Markersize',10, 'MarkerFacecolor', 'k')
hold on
plot(xi,fun(xi),'ks', 'Markersize',10)
plot(xi,gfun(xi),'ko', 'Markersize',10)
% box on
xc=(xiT(1:end-1)+xiT(2:end))/2;
d=(xiT(2:end)-xiT(1:end-1))/2;
%[inter_par,yp]=interpolateparametarization(xi,yi,sigma0./sqrt(T),1);
[inter_par]=interpolateparametarization(xi,yi,1);
[inter_par_g]=interpolateparametarization(xi,gi,1);
%x=0.01; e=max(d.^2-(x-xc).^2);
%
for ii=1:length(xx)
    x=xx(ii);
    %fr(ii)= funr(x);
    f(ii)=fun(x);
    p(ii) = interpolate_val(x,inter_par);
    pg(ii) = interpolate_val(x,inter_par_g);
    e(ii) = max(d.^2-(x-xc).^2);
    sc(ii)= (max(p(ii),pg(ii))+f0)/e(ii);
    %sc(ii) = max(p(ii)-;
    %sn(ii)= subinterpolate_search(x,e(ii),xi,yi,K);
    %sn(ii)=subquad_search(x,e(ii),xi,yi,K);
end
for ii=1:2
    x=xiu(ii);
    pU(ii) = interpolate_val(x,inter_par);
    pgU(ii) = interpolate_val(x,inter_par_g);
    eU(ii) = mindis(x,xi);
    sd(ii)= (max(pU(ii),pgU(ii))+f0)/eU(ii);
end


[tc,indc]=min(sc);
x=xx(indc);
pm = interpolate_val(x,inter_par);
pgm = interpolate_val(x,inter_par_g);
em = mindis(x,xi);
sdm= (max(pm,pgm)+f0)/em;

% %%
figure(2); clf;
plot(xx,pg,'k-','linewidth',2)
hold on
plot(xx,p,'k--', 'linewidth',2)

plot(xx,max(p,pg)+.0003, '--.', 'linewidth',4)
%  plot(xx,e*1+min(min(p,pg)),'r-')
if type =='exploring'
    plot(xx,e+min(min(p,pg))-0.08,'k-',  'linewidth',2)
 ylim([min(min(p,pg))-0.08, max(max(p,pg)+.01)])
elseif type =='replacing'
    %     plot(xx,e-0.05,'r-',  'linewidth',2)
    plot(xx,e+min(min(p,pg))-0.05,'k-',  'linewidth',2)
   ylim([min(min(p,pg))-0.05, max(max(p,pg)+.01)])
end

plot(xiu,max(pU,pgU),'k*', 'Markersize',10, 'MarkerFacecolor', 'k')
hold on
plot(xi,yi,'ks', 'Markersize',10)
plot(xi,gi,'ko', 'Markersize',10)
set(gca,'YTick',[])
set(gca,'XTick',[])
% axis square


%% replacing step Z-Omega
figure(3); clf;
subplot(2,1,1)
plot(xx, sc,'k-','linewidth',2)
hold on
% grid on
[scm, idcm]=min(sc);
[sdmm, iddm]=min(sd);
% plot(xx(idcm), scm,'ks', 'Markersize',15, 'MarkerFacecolor', 'k')

ylim([0 100])
set(gca,'YTick',[])
set(gca,'XTick',[])
% axis square
if type=='exploring'
    plot(xx(idcm), scm,'ks', 'Markersize',10,'MarkerFacecolor', 'k')
elseif type=='replacing'
    plot(xx(idcm), scm,'ks', 'Markersize',10,'MarkerFacecolor', 'k')
end

subplot(2,1,2)
% figure(4)


if type=='replacing'
    plot(xiu,sd,'k*', 'Markersize',10,'MarkerFacecolor', 'k')
hold on
% plot(xx(indc), sdm,'ks', 'Markersize',10, 'MarkerFacecolor', 'k')
plot(xx(indc), sdm,'ks', 'Markersize',10,'MarkerFacecolor', 'k')
elseif type=='exploring'
plot(xiu,sd,'k*', 'Markersize',10,'MarkerFacecolor', 'k')
hold on
% plot(xx(indc), sdm,'ks', 'Markersize',15, 'MarkerFacecolor', 'k')
plot(xx(indc), sdm,'ks', 'Markersize',10, 'MarkerFacecolor', 'k')
% plot(xx(indc), sdm,'kp', 'Markersize',15, 'MarkerFacecolor', 'k')
end

set(gca,'YTick',[])
set(gca,'XTick',[])


%% saving figures

if savekon==1
figure(1)
figure_to_publish(strcat('./figures/truth_funs_',type,'_ZOmega'))

figure(2)
figure_to_publish(strcat('./figures/inter_err_',type,'_ZOmega'))
figure(3)
figure_to_publish(strcat('./figures/searchfuns_',type,'_ZOmega'))
end





% axis square










% axis square
% grid on
%ylim([0 100])
% % errorbar(xi,yi,sigma0./sqrt(T),'k.','linewidth',2)
% %    sd = min(yp,2*yi-yp)-L*sigma0./sqrt(T);
% %[td,indd]=min(sd);
% %hold on

% %%
% %text(xx(indc),tc-.05,'x_k','fontsize',20)
% % plot(xx(indc),tc,'k*','Markersize',10)
% %text(xi(indd),td-.05,'z_k','fontsize',20)
% %set(gca,'YTickLabel',[])
% %set(gca,'XTickLabel',[])
% % ylim([-.3 10])
% grid on
%
% %
% %
% % [t,ind]=min(sn);
% % ind
% % xx(ind)
% % x=xx(ind);
% % ii=ind;
% % [y,abc]=subquad_search(x,e(ii),xi,yi,K)
% % y_inter=abc(1)*xx.^2+ abc(2)*xx+abc(3);
% % hold on
% % plot(xx,y_inter,'r-')