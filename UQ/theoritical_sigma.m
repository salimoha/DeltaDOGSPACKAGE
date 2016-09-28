function sigma=theoritical_sigma(corr,s,sigma02)
% Theoritical Calculation from corr, and sigma0
% 
%Authors:
%Shahrouz Alimo & Pooriya Beyhaghi
%March 2016
for ii=1:length(s)
    sigma(ii)=1;
    for jj=1:s(ii)-1
        sigma(ii)=sigma(ii)+2*(1-jj/s(ii))*corr(jj);
    end
    sigma(ii)=sigma02*sigma(ii)/s(ii);
end
end