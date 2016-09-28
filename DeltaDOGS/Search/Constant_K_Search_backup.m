function [x y CS]=Constant_K_Search_backup(x,inter_par, inter_Par,xc,R2)
% This funciton finds the minimizer of the search function i a simplex 
%              minimize s^k_i(x) = p^k(x) - K e^k_i(x)
%             subject to CS^k_(l,i) = g^k_l(x) - K e^k_i(x) <= 0 
%                        A x <= b 
% input:
% x: interpolating initial point for search function 
% inter_par: interpolation parameters
% xc: circumcenter for the simplex that x is located at
% R2: square of circumradiaus for the simplex that x is located at
% output:
% x: the minimizer of search function 
% y: minimum of search fucntion at x
% cse: 
% created by: Shahrouz Alimo & Pooriya Beyhaghi
% last modification: Oct/7/2015
%
% inter_Par = inter_par.inter_par_C;
% inter_par = inter_par.inter_par_S;
global n Ain bin m ms mc
% parameters of backtracking
cc=0.01;
rho=0.9; 
% Initialize the point in the simplex
%[xc,R2]=circhyp(xiT(:,tri(index,:)), n);
iter=1; cse=2;

% Calculate the Newton direction   
% search function 
y=costSearch(x,inter_par,xc, R2);
% y=cost(x,inter_par,xc, R2);
while iter<20

g=kgradSearch(x,inter_par,xc,R2);
H=khessSearch(x,inter_par,xc,R2);
H=modichol(H,0.1,20);
H=(H+H')/2;
% constrained function
CSl =[]; Jl=[]; %Hc=[];
for l=1:ms
CSl=vertcat(CSl,costSearch(x,inter_Par{l},xc, R2));
% CSl=vertcat(CSl,lkConstraint(x,inter_Par{l},xc, R2));
Jl=vertcat(Jl, kgradSearch(x,inter_Par{l},xc, R2)');
% Jl=vertcat(Jl, lkGradConstraint(x,inter_Par{l},xc, R2)');
% Hcl=khessSearch(x,inter_Par{l},xc, R2);
% Hcl=modichol(Hcl,0.1,20);
% Hcl=(Hcl+Hcl')/2;
% Hc=vertcat(Hc, Hcl);
end
% linearized constrained for sqp
A = vertcat(Ain, Jl);
b = vertcat(bin - Ain*x, -1.*CSl);
%
options=optimoptions('quadprog','Display','none');
p=quadprog(H,g,A,b,[],[],[],[],zeros(n,1),options);
% Perform the hessian modification
if norm(p)<1e-8
            break
end
a=1;
% %Backtracking ???????????????????????????
% while 1
%     x1=x+a*p;
%         y1=costSearch(x1,inter_par,xc,R2);
%         if (y-y1)>0
%         x=x1; y=y1;
%         break
%     else
%         a=a*rho;
%         if norm(a*p)<1e-4
%             break
%         end
%     end
% end
x=x+p;
iter=iter+1;
end

CS=max(CSl);
end
% % search function 
% function [M]=costSearch(x,inter_par,xc,R2)
% global K
% constant K search function 
% M=  interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2); 
% adaptive K search function 
% M=-(R2-norm(x-xc)^2)/(interpolate_val(x,inter_par)-y0); 
% if interpolate_val(x,inter_par)<y0
%     M=-inf;
% end
% end
% % kth iteration gradient for the search function
% function [ y ] = kgradSearch(x,inter_par,xc,R2)
% global  K
% % p=interpolate_val(x,inter_par);
% gp=interpolate_grad(x,inter_par);
% ge=-2*(x-xc);
% % e=R2-(x-xc)'*(x-xc);
% % constant K gradient search function 
% y = gp - K .* ge;
% %y=gp/e-(p-y0)*ge/e^2;
% % adaptive K gradient search function 
% % y=-ge/(p-y0)+e*gp/(p-y0)^2;
% end
% % kth iteration hessian for the search function
% function [ H ] = khessSearch(x,inter_par,xc,R2)
% global K n
% % p=interpolate_val(x,inter_par);
% Hp=interpolate_hessian(x,inter_par);
% % gp=interpolate_grad(x,inter_par);
% % ge=-2*(x-xc);
% % e=R2-norm(x-xc)^2;
% % Constant K hessian for the search function 
% H = Hp + K .* 2*eye(size(x,1));
% %H=Hp/e-(gp*ge.'+ge*gp.')/e^2+(p-y0)*(2*ge*ge.'/e^3+2*eye(n)/e^2);
% % adaptive K
% % H=2*eye(n)/(p-y0)+(gp*ge.'+ge*gp.')/(p-y0)^2-e*(2*gp*gp.'/(p-y0)^3-2*Hp/(p-y0)^2); 
% end
% %% Constrained functions
% function [Mc]=lkConstraint(x,inter_Par,xc,R2)
% global K m
% % constant K constrained function 
% Mc = [];
% for l = 1:m
% Mc = vertcat(Mc, [xm ym indm]= tringulation_search_constantK_constraints(inter_par_p,inter_par_g,xi,tri); - K .* (R2-norm(x-xc)^2)); 
% end 
% end
% % kth iteration gradient for the search function
% function [ yc ] = lkGradConstraint(x,inter_Par,xc,R2)
% global  K m
% yc =[];
% for l=1:m
% gp=interpolate_grad(x,inter_Par{l});
% ge=-2*(x-xc);
% % e=R2-(x-xc)'*(x-xc);
% % constant K gradient constrained function 
% yc = vertcat(yc, gp - K .* ge);
% end
% end
% % kth iteration hessian for the search function
% function [ Hc ] = lkHessConstraint(x,inter_Par,xc,R2)
% global K m n
% Hc =[];
% for l=1:m
% Hp=interpolate_hessian(x,inter_Par{l});
% % Constant K hessian for the constrained function 
% Hc = vertcat(Hc, Hp + K .* 2*eye(n));
% end
% end