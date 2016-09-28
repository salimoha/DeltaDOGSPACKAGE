% test the channel flow turulent simulation result
clear all; close all; clc;
% Initialization
% number of correlation paramters
m=18;
% tau is the correlation paramter
tau = (1:m)./(m+1);
% 
nn=2.^(9:13);
sigma_tho=nn*0; kk=1;


% load the data
disp('loading the data ....')
load('TKE_NS.mat')


% delete the transient part
disp('Removing the transient part of the signal ....')
ind = transient_time_detector(TKE);
x=TKE(ind+1:end);

% quantify the finite time averaging error for the stationary signal
disp('Quantify the finite time averaging error for the stationary signal')
%
for ii=1:length(nn)
N=nn(ii);
if length(x)<N
    break
end
%%
% normalizing the data
xx=(x(1:N)-mean(x(1:N)))./std(x(1:N));
% 
mu2_em_N = mean(xx).^2;

AR_m = 18; AR_model=ar(xx,AR_m,'burg'); a1=AR_model.a;  a1(1)=[];


[arcorr, Gama0] = ARcorr(a1,N-1,AR_model.NoiseVariance);
sigmaAR=sqrt(theoritical_sigma(arcorr,[1:N-1],Gama0))*std(x(1:N));
sigmaAR_N(ii) = sigmaAR(end);
%
[sigma2_N,theta, moment2_model, corr_model,sigmac2] = stationary_statistical_learning_reduced(x(1:N),1);
Sigma(ii)=sqrt(sigma2_N);
Mu(ii,kk)=mean(x(1:N));
sigma_theo=sigma_tho*(kk-1)/kk +emprical_sigma(x-mean(Mu(end,:)),nn)*1/kk;
end
%%
figure(1); clf;
 hold on
semilogy(0.01*nn(1:end-1),sqrt(sigma_theo(1:end-1)),'ks-','linewidth',2)
hold on
semilogy(0.01*nn(1:end-1),Sigma(1:end-1),'r-.','linewidth',2)
semilogy(0.01*nn(1:end-1),real(sigmaAR_N(1:end-1)),'b--','linewidth',2)
semilogy(0.01*nn(1:end-1),Sigma(1:end-1),'rs','linewidth',2)
semilogy(0.01*nn(1:end-1),real(sigmaAR_N(1:end-1)),'bs','linewidth',2)
grid on
grid minor
%   set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 22);
 set(gca,   'FontSize', 18);
  set(gcf, 'Position' , [60 60 700 500], 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'off' );



% 
figure(2);clf
loglog(0.01*nn(1:end-1),sqrt(sigma_theo(1:end-1)),'ks-','linewidth',2)
hold on
loglog(0.01*nn(1:end-1),Sigma(1:end-1),'r-.','linewidth',2)
loglog(0.01*nn(1:end-1),real(sigmaAR_N(1:end-1)),'b--','linewidth',2)
loglog(0.01*nn(1:end-1),Sigma(1:end-1),'rs','linewidth',2)
loglog(0.01*nn(1:end-1),real(sigmaAR_N(1:end-1)),'bs','linewidth',2)
grid on

% xlabel('T', 'fontsize',18)
% ylabel('\sigma', 'fontsize',18)
% legend( 'emperical ','Second order estimator','MLE')


grid minor
%   set(gca, 'XGrid', 'on', 'YGrid', 'on',  'FontSize', 22);
 set(gca,   'FontSize', 18);
  set(gcf, 'Position' , [60 60 700 500], 'color' , [1 1 1], 'PaperPositionMode' , 'auto', 'InverthardCopy' , 'off' );



