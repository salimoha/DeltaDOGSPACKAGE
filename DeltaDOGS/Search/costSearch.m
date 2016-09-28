function [M]=costSearch(x,inter_par,xc,R2)
global K
% constant K search function 

% M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 

if length(K) ==1
M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 
else
    M=  interpolate_val(x,inter_par) -  (K(1) + K(2)*norm(inter_par{2})).* (R2-norm(x-xc)^2); 
end
% adaptive K search function 
% M=-(R2-norm(x-xc)^2)/(interpolate_val(x,inter_par)-y0); 
% if interpolate_val(x,inter_par)<y0
%     M=-inf;
% end
end