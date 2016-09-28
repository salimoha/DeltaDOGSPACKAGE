% test the statistical analysis on kurumoto-sivashinsky &  Turbulent Channel Flow 
clear Mu Sigma

% run('../addTopath')

nn=2.^(9:13);
sigma_tho=nn*0;
%Test_CNRKW3_dtConst_4P
for kk=1:1
% filename=strcat('KS_CNRKW3', num2str(kk), '.mat');
%filename = strcat('./UncerQuant/data/02KS_CNRKW3120160223220143.mat');


%load(filename)
%x=x(t_trans+1:end);
 load 'TKE_NS.mat'
 x=TKE(ind+1:end);
%%
m=18;
tau = (1:m)./(m+1);
%%
for ii=1:length(nn)
N=nn(ii);
if length(x)<N
    break
end
%%
xx=(x(1:N)-mean(x(1:N)))./std(x(1:N));
mu2_em_N = mean(xx).^2;
AR_m = 18;
AR_model=ar(xx,AR_m,'burg');
a1=AR_model.a;  a1(1)=[];
[arcorr, Gama0] = ARcorr(a1,N-1,AR_model.NoiseVariance);
sigmaAR=sqrt(theoritical_sigma(arcorr,[1:N-1],Gama0))*std(x(1:N));
sigmaAR_N(ii) = sigmaAR(end);
%%
if kk<2
[sigma2_N,theta, moment2_model, corr_model,sigmac2] = stationary_statistical_learning_reduced(x(1:N),1);
Sigma(ii)=sqrt(sigma2_N);
end
Mu(ii,kk)=mean(x(1:N));
sigma_theo=sigma_tho*(kk-1)/kk +emprical_sigma(x-mean(Mu(end,:)),nn)*1/kk;
end
end
% %% AR test Model
% figure(1); clf;
% plot(0.2*nn,sqrt(sigma_theo),'k-','linewidth',2)
% hold on
% plot(0.2*nn,Sigma,'r.-','linewidth',2)
% plot(0.2*nn,real(sigmaAR_N),'b--','linewidth',2)
% grid on
% xlabel('N', 'fontsize',18)
% ylabel('\sigma', 'fontsize',18)

%% NS & Kuromoto Test Model
figure(2); clf;
% loglog(0.2*nn,sqrt(sigma_theo),'ks-','linewidth',2)
% hold on
% loglog(0.2*nn,Sigma,'k*--','linewidth',2)

loglog(0.2*nn(1:end-1),sqrt(sigma_theo(1:end-1)),'ks-','linewidth',2)
hold on
loglog(0.2*nn(1:end-1),Sigma(1:end-1),'r-.','linewidth',2)
loglog(0.2*nn(1:end-1),real(sigmaAR_N(1:end-1)),'b--','linewidth',2)
loglog(0.2*nn(1:end-1),Sigma(1:end-1),'rs','linewidth',2)
loglog(0.2*nn(1:end-1),real(sigmaAR_N(1:end-1)),'bs','linewidth',2)
grid on
xlabel('T', 'fontsize',18)
ylabel('\sigma', 'fontsize',18)
legend( 'emperical ','Second order estimator','MLE')