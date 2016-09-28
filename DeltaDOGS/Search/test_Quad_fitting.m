% Quad Fitting Idea
clear all
close all
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
%  xi=[min(xx) max(xx) -1 0 0.75  3 -3];
%   xi=[min(xx) max(xx) -0.925 0.75  3 -3];
% % %     xi=[min(xx) max(xx) -0.9 0.75  3 -3];
xi=[min(xx) max(xx) 0.745 0.75  3 -3];
%xi=[min(xx) max(xx) -1 0.75  3 -3];
for ii=1:length(xi)
yi(ii)=fun(xi(ii));
end
inter_par= interpolateparametarization(xi,yi,1);
for ii=1:length(xx)
yp(ii)=interpolate_val(xx(ii),inter_par);
end
figure(1); clf;
hold on
plot(xx,yy)
plot(xx,yp,'k--')
plot(xi,yi,'rs')
xlabel('x')
ylabel('f(x)')
%%
ymax = .5*max(yi); 
%ymax=200;
[H, f, d, Par] = QuadraticFitting(xi, yi, ymax, 1e-2);
a0 = [H;f;d];
q0 =0.5*a0'*Par.H*a0+ Par.f'*a0

funquad=@(x,a,b,c)(1*a*x.^2 + b*x + c);
for ii=1:length(xx)
yq(ii)=funquad(xx(ii),H,f,d);
end
plot(xx,yq,'g-')

legend('function value', 'polyharmonic spline', 'data points', 'new model')
figure(4);clf;
hold on
plot(xx,yy)
plot(xx,yp,'k--')
plot(xi,yi,'rs')
plot(xx,yq,'g-')
xlim([-5,5])
ylim([-150,200])
xlabel('x')
ylabel('f(x)')
legend('function value', 'polyharmonic spline', 'data points', 'new model')
%%


% %%
% [H, f, d, Par1] = QuadraticFitting(xi, yi, ymax, 1e2);
% a1 = [H;f;d];
% q1 =0.5*a1'*Par1.H*a1+ Par1.f'*a1
% 
% funquad=@(x,a,b,c)(1/2*a*x.^2 + b*x + c);
% for ii=1:length(xx)
% yq(ii)=funquad(xx(ii),H,f,d);
% end 
% plot(xx,yq,'g-')
% 
% 
% figure(4);clf;
% hold on
% plot(xx,yy)
% plot(xx,yp,'k--')
% plot(xi,yi,'rs')
% plot(xx,yq,'g-')
% xlim([-5,5])
% ylim([-150,200])