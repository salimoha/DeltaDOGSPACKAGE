% nn = 2.^(7:14);
nn = 2.^(20);
for i = 1:numel(nn)
%% generete the synthetic data from arima model
sigma02 = 1;mu02 = 0; A =1;
tau = [0.5, 0.25];
% tau
%use 2 AR correlation model
if length(tau) == 1
    theory_model=arima('Constant',mu02,'AR',{tau},'Variance',sigma02);
else
% % theory_model=arima('Constant',mu02,'AR',{tau(1), tau(2)},'Variance',sigma02*.1);
theory_model=arima('Constant',mu02,'AR',{tau(1), tau(2)},'Variance',sigma02*.1);
end 
N = nn(i);
s1=1:N-1;
% theory_model=arima('Constant',mu02,'AR',{-tau(1), -tau(2)},'Variance',sigma02);
x = simulate(theory_model,N); 
%
% sigma_s^2 theoretical 
tau2=theory_model.AR; tau2=-cell2mat(tau2);
% finding the correlation 
[sigma2_real, corr_real]  = theoritical_AR_sigma(tau2 , theory_model.Variance,mu02,s1);
AR_data.x = x;
AR_data.trans = 0;
AR_data.sigma2 = sigma2_real;
AR_data.corr = corr_real;
AR_data.model=theory_model;
filename = strcat('./data/Data_synt_AR',num2str(numel(tau),'%02d'),'_N_2_', num2str(log(N)/log(2),'%02d'),'.mat');  
% save(filename, 'AR_data')
end