clear all
% close all
clc
% Check the Consant Hessian Idea for Delta-Dogs:
%fun=@(x)x^4-16*x^2+5*x;
% y=-5;
x=-3;
A=40;
fun=@(y) (1.5-x+x*y)^2+(2.25-x+x*y^2)^2+(2.625-x+x*y^3)^2+A;
xx=-5:0.01:5;
% fun=@(y) y^2+A;
% fun2=@(x,y) (1.5-x+x*y)^2+(2.25-x+x*y^2)^2+(2.625-x+x*y^3)^2;
% fun=@(x)fun2
% x=1;
% fun =@(y) 100*(x - y^2)^2 + (1-y)^2;
% xx=-1.5:0.01:1.5;
for ii=1:length(xx)
yy(ii)=fun(xx(ii));
end
% xi=[min(xx) max(xx) 0.74 0.75 0.71  0 ];
xi=[min(xx) max(xx) 0.74 0.75   -3];
for ii=1:length(xi)
yi(ii)=fun(xi(ii));
end
inter_par= interpolateparametarization(xi,yi,1);
for ii=1:length(xx)
yp(ii)=interpolate_val(xx(ii),inter_par);
end
figure(1); clf;
hold on
plot(xx,yp,'k--')
plot(xi,yi,'rs')
xlabel('x')
ylabel('f(x)')

%  ylim([0 100])
%%
yiT=log(yi);
inter_par= interpolateparametarization(xi,yiT,1);
for ii=1:length(xx)
ypl(ii)=interpolate_val(xx(ii),inter_par);
end
plot(xx,yy)
figure(2); clf;
plot(xx,log(1+yy/30))
hold on
plot(xx,ypl,'k--')
plot(xi,yiT,'rs')
%%
xlabel('x')
ylabel('log(a +  b  f(x))')
