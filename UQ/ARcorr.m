function [Gama, Gama0]= ARcorr(P,m,var_eps)
% Calculate the corrolation function for the ARFIM process for lag number
% 0:n
% Calculate the covariance of ARIMA(p,d,q) 
%keyboard
[r p k]=residuez(1,[1 P]);
for i=0:m
 beta(i+1)=sum(r.*(p.^i));
end
% keyboard
Gama = xcorr(beta); Gama=Gama(m+1:2*m+1);
Gama=Gama*var_eps;
Gama0 = Gama(1);
Gama = Gama(2:end)/Gama0;
%Gama=myfilter_mex(Gama,gamma1,beta,m,mb);

