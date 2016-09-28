% 2d Islands test problem
clear all
close all
clc
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=1; m=2*n; b =1/2;
iter_max=35;
RES = 1e-2;
x_star=[0.7;0.1];
x0=[0;0];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
Method = 2;
Var = 0.2;
% end
%
% if nargin <3
% specify the objective function
%fun =@(x)x(2,:);
fun=@(x) (x(1,:)-x0(1)).^2+(x(2,:)-x0(2)).^2;
% specify the contraints functions
% con{1}=@(x) (x(2,:)-b*sin(x(1,:)*pi));
% con{2}=@(x) (-x(2,:)+b*sin(x(1,:)*pi));
%con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
%con{1}=@(x) (x(2,:)-con{1}=@(x) (x(2,:)-fun_weierstrass(x(1,:)));
%con{2}=@(x) (-x(2,:)+fun_weierstrass(x(1,:)));
KC2 = 2;
%con{1}=@(x)KC2*(x(2,:)-rastriginn(x(1,:)));
%con{2}= @(x) -con{1}(x);
con{1}=@(x) rastriginn2(x)-0.5;
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
        for k=1:iter_max
            tri=delaunayn(xi.');
      inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%%%%%%%%%%%%% FIGURE

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
   tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.6:0.6)
%     end
    plot(xi(1,:),xi(2,:),'bs') 
    if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant) ))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
%      xlim([-2 2])
%     ylim([-2 2])
    drawnow   
   
%     impro=norm(x_prev-xm);
      impro=mindis(xm,xi(:,1:end-1));

%tt = 0:0.01:1;
%plot_figures_3(xi,yi,C, tt)
%%% interpolation plot at 0
  %xx=0:0.05:1;
%   for ii=1:length(xx)
%       yg1(ii)=interpolate_val([xx(ii);0.1],inter_par_g{1});
%       yg2(ii)=interpolate_val([xx(ii);0.1],inter_par_g{2});
%       yc1(ii)=con{2}([xx(ii);0.1]);
%   end
% figure(3)
% plot(xx,yg2,'--',xx,yc1,'-')

% figure(4)
% plot(1:k,y_star,'-',1:k,-1./ym,'--')
% drawnow

%x_prev = xm;
%y_prev = ym;
% keyboard
        end