function ind = transient_time_detector(x)
% transient_time_detector(x) is an automatic procedure to determine the
% nonstationary part a signal from the stationary part.

% It finds the transient time of the simulation using the minimum variance
% intreval.
% INPUT:
% x: is the signal which after some transient part the signal becomes stationary
% OUTPUT:
% ind: is the index of signal that after that the signal could be
% considered as a stationry signal.

%Authors:
%Shahrouz Alimo & Pooriya Beyhaghi
N = length(x);
for kk = 1:floor(length(x)/2)
    y(kk) = var(x(kk+1:end))./(N-kk);
end
[~,ind] = min(y);
% keyboard
end