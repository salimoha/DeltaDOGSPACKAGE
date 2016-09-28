% parralet constraint optimization Delta-Dogs
% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints

% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K% 2d Islands test problem
% this is for constant K algorithm
clear all
close all
clc
global n m ms bnd1 bnd2 Ain bin acon Search r
r=2;
Ex=3;
% parameters
n=2; ms=1; m=2*n; b =1/2;
% number of processors
np=8;
% iter_max=26;
% % iter_max=36;
iter_max=4;
RES = 1e-2;
% x_star=[0.7;0.1]; % example 2
x_star=[0.1559 ;0.1525]; % example island
x0=[0;0];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%
% Method = 2;
% Var = 0.2;
%  Var = 0.8;
% Var = 0.02;
% % Var = 0.0476; % opt
Method = 1; %TODO
% Var =15;
% Var = 7;
 Var =12;
% Var =5;
% Var = 7;
%%
% Var = 0.3;
% KC2 = 10;
KC2 = 1;
KCf = 1;
% Method = 1;
% % Var = 0.2;
% Var = 1;
% end
%
% if nargin <3
% specify the objective function
%fun =@(x)x(2,:);
fun=@(x) ( (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2 )*KCf;
% specify the contraints functions
% con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
% con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%con{1}=@(x) (x(2,:)-con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));

%con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
%con{2}= @(x) -con{1}(x);
% con{1}=@(x) rastriginn2(x)-0.5;
% con{1}=@(x)fun_constraint;
con{1}=@(x) (rastriginn2(x)-0.5 )*KC2; %the same
%fun=@fit8;
%fun=@plygen;
%fun=@rastriginn;
% end

% if nargin <5
    lob=zeros(n,1); upb=ones(n,1);
%     lob=zeros(n,1); upb= [ 1.7;2];
% end
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
xie=xi; xiT=xi;
% iter=1
N=length(yi);
        for k=1:iter_max
            % modify the interpolations
                  inter_par_p= interpolateparametarization(xi,yi,inter_method);
             for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            % np steps of Search
            for ii=1:np
            tri=delaunayn(xi.');
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            end
            
            N=length(yi);
            
            % function evaluations
            for ii=1:np 
             yi=[yi fun(xi(:,N+ii))];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xi(:,N+ii))];
            end
            end

   h=figure(2); %clf;
%    tri=delaunayn(xi.');
%     triplot(tri,xi(1,:),xi(2,:))
    hold on
   tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.5:0.5)
colormap gray
colormap bone
brighten(0.5)
%     end
strmax = ['   ', num2str(k)];
% text(xi(1,:),xi(2,:),strmax ,'HorizontalAlignment','right');
    plot(xi(1,end),xi(2,end),'rs', 'MarkerSize', 15) 
      plot(xi(1,:),xi(2,:),'rs', 'MarkerSize', 15) 
    if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant) ))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
%      xlim([-2 2])
%     ylim([-2 2])
%    triplot(tri,xi(1,:),xi(2,:))
drawnow   
% % % 
% % % %     impro=norm(x_prev-xm);
% % %       impro=mindis(xm,xi(:,1:end-1));
% % % %   saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
tt = 0:0.01:1;

%%
impro=mindis(xm,xi(:,1:end-1));
if impro<RES|| k == iter_max+1
    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(33));
    % keyboard
%      save(strcat(fName,'_workspace','.mat')) 
%   csvwrite(strcat(fName,'.csv'), iter_1);
   break
end






        end
        

        