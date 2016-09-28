function [x y CS]=Constant_K_Search_interior(x,inter_par, inter_Par,xc,R2,Search)
% This funciton finds the minimizer of the search function i a simplex 
%              minimize s^k_i(x) = p^k(x) - K e^k_i(x)
%             subject to CS^k_(l,i) = g^k_l(x) - K e^k_i(x) <= 0 
%                        A x <= b 
% input:
% method: interior point
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

global n Ain bin ms  K Eps
Eps = 2e-3;
K = Search.constant;
% parameters of backtracking
rho=0.9; 
% Initialize the point in the simplex
%[xc,R2]=circhyp(xiT(:,tri(index,:)), n); 

% Calculate the Newton direction   
% search function 

% Calculate an initail feasible point
[x,CS]=Constrained_feasible_point_finder(x,inter_Par,xc,R2);
% keyboard
if CS>Eps
    y=inf;
else
% y=cost(x,inter_par,xc, R2);
% mu=0.01;
mu=1e-6;
for iterrr=1:10
iter=1;
% solve for given mu
%y=0.01*ones(ms,1);
for jj=1:ms
    c(jj)=-interpolate_val(x,inter_Par{jj})+K .* (R2-norm(x-xc)^2);
end
keyboard
% y=mu./(c)';
x_pre=x;
%keyboard
while iter<3
% Calculate the primial_dual_search function and its derivatives 
M=interior_cost(x,inter_par,inter_Par,mu,xc,R2);
DM=interior_grad(x,inter_par,inter_Par,mu,xc,R2);
D2M=interior_hessian(x,inter_par,inter_Par,mu,xc,R2);
%D2M=primial_dual_hessian(x,y,inter_par,inter_Par,xc,R2);
D2M=modichol(D2M,0.1,20);
D2M=(D2M+D2M')/2;

% quadprog
options=optimoptions('quadprog','Display','none');
% pp=quadprog(D2M,DM,[Ain zeros(2*n,ms); zeros(ms,n) -eye(ms)],[bin - Ain*x; y],[],[],[],[],zeros(n+ms,1),options);
 p=quadprog(D2M,DM,Ain,[bin - Ain*x],[],[],[],[],zeros(n+ms,1),options);

% backtracking
a=1;
while 1
    x1=x+a*p;
        M1=interior_cost(x1,inter_par,inter_Par,mu,xc,R2);
        if (M-M1)>0
        x=x1; M=M1;
             break
        else
        a=a*rho;
           if norm(a*[p ;q])<1e-4
             break
           end
        end
end
if norm(a*p)<1e-12
    break
end

iter=iter+1;
end
impro=norm(x_pre-x);
if impro<1e-4
    break
end
mu=mu/2;
end

CS=-Inf;
for l=1:ms
CS=max([CS,interpolate_val(x,inter_Par{l}) - K .* (R2-norm(x-xc)^2)]);
end
y=interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2);
end

end

function [M]=interior_cost(x,inter_par,inter_Par,mu,xc,R2)
global ms K 
M=interpolate_val(x,inter_par) - K .* (R2-norm(x-xc)^2);
for jj=1:ms
    c=max(0,-interpolate_val(x,inter_Par{jj})+K .* (R2-norm(x-xc)^2));
%     y(jj)=max(0,y(jj));
    M=M-mu*log(c);
end
end


function [DM]=interior_grad(x,inter_par,inter_Par,mu,xc,R2)

global ms K
DM=interpolate_grad(x,inter_par)+2*K*(x-xc);
for jj=1:ms
    c=-interpolate_val(x,inter_Par{jj})+K*(R2-norm(x-xc)^2);
    J=-interpolate_grad(x,inter_Par{jj})-2*K*(x-xc);
    DM=DM-mu*J/c;
end
%DM=[g-J'*(2*pii-y);c-mu./y];
end

%function [D2M]=primial_dual_hessian(x,y,inter_par,inter_Par,R2,mu)
function [D2M]=interior_hessian(x,inter_par,inter_Par,mu, xc,R2)
global ms K
Hm=interpolate_hessian(x,inter_par)+2*K*eye(length(x));
for jj=1:ms
    c=-interpolate_val(x,inter_Par{jj})+K*(R2-norm(x-xc)^2);
    J=-interpolate_grad(x,inter_Par{jj})+2*K*(x-xc);
    H=-interpolate_hessian(x,inter_Par{jj})-2*K*eye(length(x));
    Hm = Hm - mu/c*H + mu / c^2 * J * J';
    %      Hm=Hm-2*y(jj)*H;
%     Hm=Hm-y(jj)*H;
end
D2M = Hm;
%  D2M=[Hm+2*J'*diag(y./(c)')*J J'; J diag((c)'./y)];
 %D2M=[Hm+2*J'*diag(mu./(c.^2)')*J J'; J diag(mu./y.^2)];
end
