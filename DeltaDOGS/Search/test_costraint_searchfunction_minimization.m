close all
clc
index = 3
[xc,R2]=circhyp(xi(:,tri(index,:)), n); 
xv=0:0.1:1;
K=1;
for ii=1:length(xv)
    for jj=1:length(xv)
        x=[xv(ii); xv(jj)];
        e=R2-norm(x-xc)^2;
        S(jj,ii)=interpolate_val(x,inter_par_p)-K*e;
        g1(jj,ii)=interpolate_val(x,inter_par_g{1})-K*e;
        g2(jj,ii)=interpolate_val(x,inter_par_g{2})-K*e;
    end
end
%%
figure(1)
contourf(S)
% colorbar
figure(2)
 contourf(xv,xv,max(g1,g2),-1:1:0)

% contourf(xv,xv,g1,-1:1:0)
% legend(strcat(' K = ',  , 

%%%%%%%%%%%%


xv=0:0.1:1;
K=2;
for ii=1:length(xv)
    for jj=1:length(xv)
        x=[xv(ii); xv(jj)];
        e=R2-norm(x-xc)^2;
        S(jj,ii)=interpolate_val(x,inter_par_p)-K*e;
        g1(jj,ii)=interpolate_val(x,inter_par_g{1})-K*e;
        g2(jj,ii)=interpolate_val(x,inter_par_g{2})-K*e;
    end
end

figure(3)
contourf(S)
% colorbar
figure(4)
contourf(xv,xv,max(g1,g2),-1:1:0)


