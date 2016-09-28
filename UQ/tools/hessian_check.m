function [ hessian_numeric, hessian_theory, err] = hessian_check(x,fun,delta,hess)
% calculates the finite difference gradiant to check the gradient formula
n = length(x);
% keyboard
y0 = fun(x);
          if nargin < 4
              err = nan;
              hessian_theory = nan(n,n);
          end             
for ii=1:n
   for jj = 1:n
        Delta_i = [zeros(ii-1,1) ; 1; zeros(n-ii,1)];
         Delta_j = [zeros(jj-1,1) ; 1; zeros(n-jj,1)];
           if length(y0)>1; 
      hessian_numeric(ii,jj) = sum((fun(x+Delta_i.*delta +Delta_j.*delta ) - ...
            fun(x+Delta_j.*delta )-fun(x+Delta_i.*delta) + fun(x))./delta);
           else
                 hessian_numeric(ii,jj) =((fun(x+Delta_i.*delta +Delta_j.*delta ) - ...
            fun(x+Delta_j.*delta )-fun(x+Delta_i.*delta) + fun(x))./delta^2);
               
           end
   end
   
end

     if nargin == 4
               hessian_theory = hess(x);
                err = norm(hessian_numeric - hessian_theory);
     end

end
