function corr=exponential_correlation(A,tau,N)
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