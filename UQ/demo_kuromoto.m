% test the statistical analysis on kurumoto-sivashinsky
clear Mu Sigma
nn=2.^(8:14);
sigma_tho=nn*0;
%Test_CNRKW3_dtConst_4P
for kk=1:1
filename = strcat('./KS_SIM_NEW_1024');
load(filename)
x=TKE;
x=x(t_trans+1:end);
for ii=1:length(nn)
N=nn(ii);
if length(x)<N
    break
end
if kk<2
[sigma2_N,theta, moment2_model, corr_model,sigmac2] = stationary_statistical_learning_reduced(x(1:N),1);
Sigma(ii)=sqrt(sigma2_N);
end
Mu(ii,kk)=mean(x(1:N));
sigma_theo=sigma_tho*(kk-1)/kk +emprical_sigma(x-mean(Mu(end,:)),nn)*1/kk;
end
end
%%
figure(1); clf;
plot(0.2*nn,sqrt(sigma_theo),'ks-','linewidth',2)
hold on
plot(0.2*nn,Sigma,'k*--','linewidth',2)
grid on
xlabel('N', 'fontsize',24)
ylabel('\sigma', 'fontsize',24)

figure(2); clf;
% loglog(0.2*nn,sqrt(sigma_theo),'ks-','linewidth',2)
% hold on
% loglog(0.2*nn,Sigma,'k*--','linewidth',2)

loglog(0.2*nn(1:end-1),sqrt(sigma_theo(1:end-1)),'ks-','linewidth',2)
hold on
loglog(0.2*nn(1:end-1),Sigma(1:end-1),'k*--','linewidth',2)
grid on
xlabel('N', 'fontsize',16)
ylabel('\sigma^2', 'fontsize',16)
legend( 'emperical ','new method', 'location', 'best')
% set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 22);
 set(gca,   'FontSize', 18);
% set(gcf, 'Position' , [60 60 700 500], 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'off' );
% set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
% set(gcf, 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'on' );
%  set(gcf,'units','normalized', 'outerposition', [0,0,0.5,0.5])

figure(3); clf;
% nn(end)=[]; Sigma(end)=[]; Mu(end)=[];
semilogx(nn, mean(Mu(end,:))+ Sigma,'k--','linewidth',2)
hold on
semilogx(nn, Mu,'s-','linewidth',2)
semilogx(nn, mean(Mu(end,:))- Sigma,'k--','linewidth',2)
 set(gca,   'FontSize', 18);
grid on
xlabel('N', 'fontsize',24)
ylabel('\mu', 'fontsize',24)
%mu_ave=mean(Mu(end,:));
figure(4); clf;
plot(nn, mean(Mu(end,:))+ Sigma,'k--','linewidth',2)
hold on
plot(nn, Mu,'s-','linewidth',2)
plot(nn, mean(Mu(end,:))- Sigma,'k--','linewidth',2)
grid on
xlabel('N', 'fontsize',24)
ylabel('\mu', 'fontsize',24)
%mu_ave=mean(Mu(end,:));


% 
% 
% close all
% figure(1)
% loglog(0.2*nn,sqrt(sigma_theo),'k-','linewidth',2)
% hold on
% loglog(0.2*nn,Sigma,'k--','linewidth',2)
% 
% figure(2)
% semilogx(nn, mean(Mu(end,:))+ Sigma,'k--')
% hold on
% semilogx(nn, Mu,'-')
% hold on
% semilogx(nn, mean(Mu(end,:))- Sigma,'k--')
% %mu_ave=mean(Mu(end,:));
% 
% 
% 
% 
