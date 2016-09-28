clear all; close all; clc;
warning('off','all')
% build orginal theoretical data
%  keyboard
sigma02 = 1;mu02 = 0;
%A =[0.5;0.5];
%tau = [0.5;0.25];
% different AR models
p=poly([0.1 0.95 0.8 0.7 0.5*exp(1i*pi*0.3) 0.5*exp(-1i*pi*0.3)]); p(1)=[];
% Take the first five values = C, where C equals
% (a) 0
% (b) 100
CC=[100];
n=9; N=2^n; 
Simulate_Length=1000;
ind_mean=0;
ind_COV=0;
DATA_mean=zeros(N+1,1);
DATA_COV=zeros(N+1,1);
IND=[];
    for jj=1:length(CC)
        C=CC(jj);
        for kk=1:Simulate_Length
            x0 = ones(length(p),1)*C;
            x = ARsimulator(p,N,x0,sigma02*0.0);
            x = ARsimulator(p,N,x0,sigma02*0.1);
            ind = transient_time_detector(x);
            IND=[IND ind];
           % keyboard
           % ind_mean1=ind_mean*(kk-1)/kk+ind/kk;
            DATA_mean1 = DATA_mean*(kk-1)/kk+ x*1/kk;
            if kk>1
            
            DATA_COV = ((kk-2)*DATA_COV + (kk-1)*(DATA_mean1-DATA_mean).^2+(x-DATA_mean).^2)/(kk-1);
           % ind_COV = ((kk-2)*ind_COV + (kk-1)*(ind_mean1-ind_mean)^2+(ind-ind_mean)^2)/(kk-1);
            end
            DATA_mean=DATA_mean1;
        %    ind_mean=ind_mean1;
        end
%         DATA_low=DATA_mean+sqrt(DATA_COV);
%         DATA_up=DATA_mean+sqrt(DATA_COV);
%         ind_low=ind_mean+sqrt(ind_COV);
%         ind_up=ind_mean+sqrt(ind_COV);
        figure(1)
            plot(DATA_mean,'k-', 'linewidth',2)
            hold on
            plot(DATA_mean+sqrt(DATA_COV),'k--', 'linewidth',1.5)
            plot(DATA_mean-sqrt(DATA_COV),'k--', 'linewidth',1.5)
            plot(DATA_mean,'k-', 'linewidth',2)

            set(gca,'fontsize', 18)
            xlim([0 500])
%             print Mean_SIG_100_AR.eps
        figure(2)
            [f,x]=hist(IND,50);
            bar(x,f/sum(f));
            xlim([0 500])
%             print -depsc HIS_0_AR.eps
            %U= [-60 60];
            %plot([ind_mean ind_mean], U,'k-')
            %plot([ind_mean+sqrt(ind_COV) ind_mean+sqrt(ind_COV)], U,'k--')
            %plot([ind_mean-sqrt(ind_COV) ind_mean-sqrt(ind_COV)], U,'k--')
            
        
    end
    %%
    
    
    close all
%     ciplot(,DATA_mean+sqrt(DATA_COV))
lowerPI=DATA_mean-sqrt(DATA_COV)*1.;
upperPI=DATA_mean+sqrt(DATA_COV)*1.;
t=1:length(upperPI);
qq =[ 1:10:90 ]/100;
i = floor(length(qq));

kk=1;
cmap = [0.9 0.9*i/length(qq) 0.9*i/length(qq)];
ciplot(lowerPI(kk:end), upperPI(kk:end),t(kk:end), cmap);
%        end
hold on,

x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
plot(x,'-.', 'linewidth',0.9)
x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
plot(x,'b:', 'linewidth',0.9)
x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
plot(x,'--', 'linewidth',0.6)
x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
plot(x,'-.', 'linewidth',0.6)
plot(DATA_mean(kk:end),'k-', 'linewidth',2)
% %%%%%
% x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
% plot(x,'-.', 'linewidth',0.9)
% x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
% plot(x,':', 'linewidth',0.9)
% x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
% plot(x,'--', 'linewidth',0.6)
% x = ARsimulator(p,N,x0,sigma02*0.1);x=x(kk:end);
% plot(x,'-.', 'linewidth',0.6)
% %%%%%

xlim([7 500])
ylim([-90 120])
figure_to_publish('./Mean_SIG_100_AR')
% ylim([-80,80])
% figure_to_publish('./Mean_SIG_0_AR')
% grid minor

% figure_to_publish('./Mean_SIG_100_AR')