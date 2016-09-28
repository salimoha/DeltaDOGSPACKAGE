function [sigma]=emprical_sigma(x,s)
% emprical Calculation of sigma(s).
% sigma= \sum mean(x(t:t+s)).^2
% input:
% x: data
% s: desired averaging lengths
% output: sigma
%  keyboard
% %%
% mu=[];
N=length(x);
for jj=1:length(s)
for i=0:floor(N/s(jj))-1
   mu(i+1)=mean(x(i*s(jj)+1:(i+1)*s(jj))); 
%    mu=[mu, mean(x(i*s(jj)+1:(i+1)*s(jj)))]; 
end
sigma(jj)=mean(mu.^2);
clear mu
% mu=[];
end
% keyboard
end