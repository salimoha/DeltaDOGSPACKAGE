function [k] = consraints_delta_dogs(Var, Method, fun, con, lob, upb)
% The Delta_Dogs algorithm with nonlinear costraints.
% Simple bound costraints for the control parameter is imposed.
% bnd1<x<bnd2

% fun: The objective function
% con{i}= The costraint function

% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints

% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K
% clear all
% close all
% clc
global n m ms bnd1 bnd2 Ain bin acon Search r 
r=5;
% parameters
n=2; ms=2; m=2*n;
if nargin <2 
% Search.method=2;
% Search.constant=0;
Method = 2;
Var = 0;
end
%
if nargin <3
% specify the objective function
fun =@(x) 1-x(2,:);
% specify the contraints functions
con{1}=@(x) x(2,:)-sin(x(1,:)*pi);
con{2}=@(x) -x(2,:)+sin(x(1,:)*pi);
%fun=@fit8;
%fun=@plygen;
%fun=@rastriginn;
end

if nargin <4
    bnd1=zeros(n,1); bnd2=ones(n,1); 
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
%ss0=[0.01 1.0472 0.01 0.5 0.1 0.7 0.3 0 0 0].';
%ss1=[pi/2 pi/2 pi/2 1.5 0.5 1.5 0.5 1 1 1].';
%ss0=[10 1e-4 10 1e-4 10 1e-4 100 1e-4]';
%ss1=[1000 1 100 2.5 20 1.6 10000 pi/2]';
% the inequality constraints
%Ain=[]; bin=[];
% the equlaity constraints
%Aeq=[]; beq=[];
% initial calculations
% bnd1=zeros(n,1); bnd2=ones(n,1); 
%rk=ss1-ss0;
% input the inequality constraints
%Ain=[]; bin=[];
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
  xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
x_prev = xi(:,1);
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
delta_tol=0.2;
iter_max=50;
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

  figure(1); clf;
    subplot(2,1,1)
    plot(1:length(yi(1:end)),yi,'-')
    title('objective function')
    grid on
    subplot(2,1,2)
    plot(1:length(yi(1:end)),C{1},'-')
    title('costrained function')
    grid on
   figure(2); clf;
   tri=delaunayn(xi.');
    triplot(tri,xi(1,:),xi(2,:))
    hold on
%     tt = 0:0.01:1;
    
    tt = -2:0.01:2;
    global EX KC2
    %%
    switch EX
    case 1
%     ex1
%     plot(tt, 1*sin(pi*tt), 'r', 'linewidth', 3)
      plot(tt, fun_weierstrass(tt), 'r', 'linewidth', 3)
%    ex2        
        case 2
      plot(tt, 1/2*fun_weierstrass(tt), 'r', 'linewidth', 3)
    
        case 3
         plot(tt, KC2*fun_weierstrass(tt), 'r', 'linewidth', 3)     
        case 4
         plot(tt, KC2*fun_weierstrass(tt), 'r', 'linewidth', 3)     
            
    end
    plot(xi(1,:),xi(2,:),'s')
    
    if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant) ))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
%     xlim([0 1])
%     ylim([0 1])
     xlim([-1.7 1.7])
    ylim([-2 2])
    drawnow   
    
%     keyboard
%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));
if impro<1e-4 || k == iter_max

    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(EX));
%      keyboard
     save(strcat(fName,'_workspace','.mat')) 
  csvwrite(strcat(fName,'.csv'), iter_1);
    break
end
x_prev = xm;
        end