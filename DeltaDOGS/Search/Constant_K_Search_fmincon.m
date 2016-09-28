function [x y CS]=Constant_K_Search_fmincon(x,inter_par_p, inter_par_g,xc,R2,Search)
% This funciton finds the minimizer of the search function i a simplex 
%              minimize s^k_i(x) = p^k(x) - K e^k_i(x)
%             subject to CS^k_(l,i) = g^k_l(x) - K e^k_i(x) <= 0 
%                        A x <= b 
% input:
% x: interpolating initial point for search function 
% inter_par_p: interpolation parameters
% xc: circumcenter for the simplex that x is located at
% R2: square of circumradiaus for the simplex that x is located at
% output:
% x: the minimizer of search function 
% y: minimum of search fucntion at x
% cse: 
% created by: Shahrouz Alimo & Pooriya Beyhaghi
% last modification: Oct/7/2015
%
global n Ain bin ms  K
K = Search.constant;
% parameters of backtracking
rho=0.9; 
% Initialize the point in the simplex
%[xc,R2]=circhyp(xiT(:,tri(index,:)), n); 

% Calculate the Newton direction   
% search function 

% Calculate an initail feasible point
[x,CS]=Constrained_feasible_point_finder(x,inter_par_g,xc,R2);
if CS>0
    y=inf;
else
% y=cost(x,inter_par_p,xc, R2);
mu=1;
for iterrr=1:15
iter=1;
% solve for given mu
y=0.01*ones(ms,1);
x_pre=x;
while iter<10
% Calculate the primial_dual_search function and its derivatives 
M=primial_dual_cost(x,y,inter_par_p,inter_par_g,mu,xc,R2);
DM=primial_dual_grad(x,y,inter_par_p,inter_par_g,mu,xc,R2);
D2M=primial_dual_hessian(x,y,inter_par_p,inter_par_g,xc,R2);
D2M=modichol(D2M,0.1,20);
D2M=(D2M+D2M')/2;

% quadprog
options=optimoptions('quadprog','Display','none');
pp=quadprog(D2M,DM,[Ain zeros(2*n,ms); zeros(ms,n) -eye(ms)],[bin - Ain*x; -y],[],[],[],[],zeros(n+ms,1),options);
p=pp(1:n); q=pp(n+1:n+ms);
% backtracking
a=1;
while 1
    x1=x+a*p;
    y1=y+a*q;
        M1=primial_dual_cost(x1,y1,inter_par_p,inter_par_g,mu,xc,R2);
        if (M-M1)>0
        x=x1; y=y1; M=M1;
             break
        else
        a=a*rho;
           if norm(a*[p ;q])<1e-4
             break
           end
        end
end
if norm(a*p)<1e-4
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
CS=max([CS,interpolate_val(x,inter_par_g{l}) - K .* (R2-norm(x-xc)^2)]);
end
y=interpolate_val(x,inter_par_p) - K .* (R2-norm(x-xc)^2);
end

%%

options = optimset('Algorithm','interior-point',...
'Display','iter','GradObj','on','GradConstr','on',...
'Hessian','user-supplied','HessFcn',@hessinterior);

[x fval mflag output]=fmincon(@ObjSearchFun,x,...
[],[],[],[],[],[],@NonliCon,options)

end


%%  fmin functions
function [M DM] = ObjSearchFun(x)
% Calculate objective f
global inter_par_p K xc R2
M = interpolate_val(x, inter_par_p) -K * (R2-norm(x-xc)^2);
if nargout > 1 % gradient required
    DM = interpolate_grad(x,inter_par_p)  + 2*K*(x-xc);

end
end

function [c J] = NonliCon(x)
% Calculate state constraints values and grads
global inter_par_g K xc R2 ms
    for jj=1:ms
    c(jj)=-interpolate_val(x,inter_par_g{jj})+K*(R2-norm(x-xc)^2);
    end
    
    if nargout >1
     for jj=1:ms
    J(jj,:)=-interpolate_grad(x,inter_par_g{jj})+2*K*(x-xc);
     end
    end
end
%% hessian

function h = hessinterior(x,lambda)
% calculates the hessian of search fucntions and lagrangian constraints
global inter_par_g K xc R2 ms
h = interpolate_hessian(x,inter_par_p)+2*K*eye(length(x));% Hessian of search function

for jj=1:ms
    hessc(jj,:,:) = interpolate_hessian(x,inter_par_g{jj})+2*K*eye(length(x));
end
for jj =1:ms
h = h + lambda.ineqnonlin(jj)*hessc(jj,:,:);
end
end


%%




function [M]=primial_dual_cost(x,y,inter_par_p,inter_par_g,mu,xc,R2)
global ms K
M=interpolate_val(x,inter_par_p) - K .* (R2-norm(x-xc)^2);
for jj=1:ms
    c(jj)=max(0,-interpolate_val(x,inter_par_g{jj})+K .* (R2-norm(x-xc)^2));
    M=M-mu*log(c(jj))-mu*(log(c(jj)*y(jj)/mu)+(1-c(jj)*y(jj)/mu));
end
end


function [DM]=primial_dual_grad(x,y,inter_par_p,inter_par_g,mu,xc,R2)
global ms K
g=interpolate_grad(x,inter_par_p)+2*K*(x-xc);
for jj=1:ms
    c(jj)=-interpolate_val(x,inter_par_g{jj})+K*(R2-norm(x-xc)^2);
    J(jj,:)=-interpolate_grad(x,inter_par_g{jj})-2*K*(x-xc);
end
c=c';
pii=mu./c;
D=c./y;
DM=[g-J'*(2*pii-y);(y-pii).*D]; 
end

function [D2M]=primial_dual_hessian(x,y,inter_par_p,inter_par_g,xc,R2)
global ms K
Hm=interpolate_hessian(x,inter_par_p)+2*K*eye(length(x));
for jj=1:ms
    c(jj)=-interpolate_val(x,inter_par_g{jj})+K*(R2-norm(x-xc)^2);
    J(jj,:)=-interpolate_grad(x,inter_par_g{jj})+2*K*(x-xc);
    H=-interpolate_hessian(x,inter_par_g{jj})-2*K*eye(length(x));
    Hm=Hm-2*y(jj)*H;
end

 D2M=[Hm+2*J'*diag(y./(c)')*J J'; J diag((c)'./y)];
end