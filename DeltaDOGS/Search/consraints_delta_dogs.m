function [k] = consraints_delta_dogs(Var, Method, fun, con,lob,upb)
% The Delta_Dogs algorithm with nonlinear costraints.
% Simple bound costraints for the control parameter is imposed.
% bnd1<x<bnd2

% fun: The objective function
% con{i}= The costraint function

% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints

% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K
if nargin <1
clear all
close all
clc
end
global n m ms bnd1 bnd2 Ain bin acon Search r 
r=2;
% parameters
n=2; ms=2; m=2*n; b =1/2;
iter_max=30;
if nargin <2 
% Search.method=2;
% Search.constant=0;
Method = 2;
Var = 0;
 end
%
 if nargin <3
% specify the objective function
fun =@(x)x(2,:);
% specify the contraints functions
con{1}=@(x)100*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
% con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
% con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%con{1}=@(x) (x(2,:)-con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
end

if nargin <5
    lob=zeros(n,1); upb=ones(n,1);
%     lob=zeros(n,1); upb= [ 1.7;2];
 end
%
Search.method = Method;
Search.constant = Var;
bnd1 = lob;
bnd2 = upb;
% Search.method=1;
% Search.constant=1;
% interpolaion strategy
inter_method=1;
% upper and lower bouds

% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
  xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
for ii=1:m
    acon{ii}=[];
end
for ii=1:size(xi,2)
    ind=1:2*n;
    ind=ind(Ain*xi(:,ii)-bin>-0.01);
    for jj=1:length(ind)
    acon{ind(jj)}=[acon{ind(jj)} ii];
    end
    yi(ii)=fun(xi(:,ii));
    for jj=1:ms
        C{jj}(ii)=con{jj}(xi(:,ii));
    end
end
x_prev = xi(:,1);
y_prev = yi(:,1);
delta_tol=0.2;

xi(1,1)=xi(1,1)+0.00001;
        for k=1:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
            % feasible constraint projection
            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end

% ploting the steps and convergance 
            tt = 0:0.01:1;
     plot_figures_3(xi,yi,C, tt)    
%%% ploting interpolation gi and constraint function ci
      xx=0:0.05:1;
  for ii=1:length(xx)
      yg1(ii)=interpolate_val([xx(ii);0],inter_par_g{1});
      yg2(ii)=interpolate_val([xx(ii);0],inter_par_g{2});
      yc1(ii)=con{2}([xx(ii);0]);
  end
   figure(3)
plot(xx,yg2,'--',xx,yc1,'-', 'linewidth', 3)  
legend('intepolation constraint' , 'analytical constraint')

%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
if impro<1e-2 || k == iter_max

    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(33));
     %keyboard
%      save(strcat(fName,'_workspace','.mat')) 
%   csvwrite(strcat(fName,'.csv'), iter_1);
    break
end


x_prev = xm;
y_prev = ym;
 keyboard
        end
        
       
        