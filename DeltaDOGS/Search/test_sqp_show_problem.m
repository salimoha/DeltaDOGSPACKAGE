xv=0:0.005:1;
clear U U_l U_0 xr jj ii
for ii=1:length(xv)
    for jj=1:length(xv)
        xr=[xv(ii) ;xv(jj)];
        U(ii,jj)=costSearch(xr,inter_Par{2},xc,R2);
        U_l(ii,jj)=CSl(2)+ Jl(2,:)*(xr-x);
        U_0(ii,jj)=0;
    end
end
figure(2)
clf
contour(xv,xv,U', [0,0],'r-')
colormap(autumn)
hold on
contour(xv,xv,U_l', [0,0],'k-')
colormap(winter)
hold on
contour(xv,xv,U_0,[0,0])
colormap(summer)

%%
figure(2)
clf
contour(xv,xv,U', [0,0])
colormap(autumn)
hold on
surf(U_l')
colormap(winter)
hold on
surf(U_0)
colormap(summer)



