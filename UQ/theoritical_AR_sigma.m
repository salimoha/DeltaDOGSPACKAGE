function  [sigma2, corr] =theoritical_AR_sigma(tau , sigma02,mu2,s,var_eps)
% Theoritical Calculation of sigma(s).
% input:
% x: data
% s: desired averaging lengths
% output: sigma
% if nargin <5
%     var_eps = 0.1;
% end
%keyboard
% corr = ARIMAcorr(tau,[],0,length(s)-1);
[corr, Gama0] = ARcorr(tau,length(s)-1,var_eps);
keyboard
sigma2=theoritical_sigma(corr,s,Gama0);
sigma2=mu2+sigma2;
end