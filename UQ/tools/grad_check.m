function [ grad_numeric] = grad_check(x,fun,delta)
% calculates the finite difference gradiant to check the gradient formula
 n = length(x);
% keyboard
err=0;
for ii=1:n
    Delta_k = [zeros(ii-1,1) ; 1; zeros(n-ii,1)];
     grad_numeric(ii,:) = (fun(x+Delta_k.*delta) - fun(x))./delta;
%    grad_numeric(ii) = (fun(x(ii)+delta) - fun(x(ii)))./delta;
%     grad_theory(ii) = grad(x(ii));
%     err = err + ( grad_theory(ii) -  grad_numeric(ii,ii))^2;
end


end

