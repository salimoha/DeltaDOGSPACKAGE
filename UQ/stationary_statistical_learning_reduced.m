function [sigma2_N,theta, moment2_model, corr_model,sigmac2] = stationary_statistical_learning_reduced(x,m)
% this fucntion gets the variables as:
% INPUT:
% x: stationary ergodic process signal 
% m: moment order for the process (number of correlation parameters)
% cost fucntion
% [L,gradL, HessL] =Loss_fun(sigmaT2,mu2, A, tau, sigmac2,N,M)
%% finding optimum theta value
% keyboard

N=length(x);
M=floor(sqrt(N));
s=1:2*M;
variance=var(x);
x=(x-mean(x))./std(x);


sigmac2=emprical_sigma(x,s);
tau=(1:m)./(m+1);

% find optimum theta
options = optimoptions('fmincon','display','none','GradObj','on');
fun=@(tau) Loss_fun_reduced(tau,sigmac2);
theta(m+1:2*m)=fmincon(fun,tau,[],[],[],[],zeros(m,1),ones(m,1),[],options);

theta(1:m)=optimum_A(theta(m+1:2*m),sigmac2);
% find modeled sigma2
% [moment2_model, corr_model] =spMoment2([1 0 theta],N);
[moment2_model, corr_model] =Thoe_moment2([1 0 theta],N);

sigma2_N=moment2_model(:,end)*variance;
end

function [L,DL] =Loss_fun_reduced(tau,sigmac2)
% Reduced wighted least square.
% tau is fixed.

L = 0;
%keyboard
m=length(tau);
H=zeros(m);
for ss = 1:length(sigmac2)
for ii=1:m
as=1:ss;
Ls(ii)=1/ss*(1+2*(1-as./ss)*tau(ii).^as')-sigmac2(ss);
DL(ii)=2*Ls(ii)*2/ss*((as-as.^2/ss)*tau(ii).^(as-1)');
end
H=H+Ls'*Ls;
end
%keyboard
options=optimoptions('quadprog','Algorithm','active-set','Display','none');
[A,L] =quadprog(H,zeros(m,1),-eye(m),zeros(m,1),ones(1,m),1,[],[],ones(m,1)/m,options);
DL=(DL'.*A);
end



function A =optimum_A(tau,sigmac2)
% Reduced wighted least square.
% tau is fixed.

L = 0;
%keyboard
m=length(tau);
H=zeros(m);
for ss = 1:length(sigmac2)
for ii=1:m
as=1:ss;
Ls(ii)=1/ss*(1+2*(1-as./ss)*tau(ii).^as')-sigmac2(ss);
DL(ii)=2*Ls(ii)*2/ss*((as-as.^2/ss)*tau(ii).^(as-1)');
end
H=H+Ls'*Ls;
end
%keyboard
options=optimoptions('quadprog','Algorithm','active-set','Display','none');
A=quadprog(H,zeros(m,1),-eye(m),zeros(m,1),ones(1,m),1,[],[],ones(m,1)/m,options);
end




function  [moment2, corr_model] =Thoe_moment2(theta,N)
% Theoritical Calculation of moment(s).
% input:
% x: data
% s: desired averaging lengths
% output: moment2
%keyboard
% corr = ARIMAcorr(tau,[],0,length(s)-1);
m = length(theta)/2-1;
corr_model= exponential_correlation( theta(3:m+2),theta(m+3:end),N);
s = [1:N]';
moment2=theta(2) + theoritical_sigma(corr_model',s,theta(1));
end


function corr= exponential_correlation(A,tau,N)
% Theoritical Calculation of sigma(s).
% input:
% x: data
% s: desired averaging lengths
% output: sigma
%keyboard
corr=zeros(1,N);
for ii=1:length(tau)
corr=corr+A(ii)*tau(ii).^(1:N);
end
end

