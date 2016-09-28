% Local optimization method based on Delaunay triangulayion
clear all
close all
%%
fun=@(x) sum(x.^2);
%fun=@(x) (x(1)-1)^2+10*(x(2)-x(1)^2)^2;
xi=2*[-1 -1 1 1 -1; -1 1 -1 1 0.5];
for ii=1:size(xi,2)
    yi(ii)=fun(xi(:,ii));
end
y0=-0.1;
%y0=-500;
K=10000;
% K=.1;
for kk=1:20
    clear ss
    figure(1); clf;
[ymin,ind_min]=min(yi);
tri=delaunayn(xi.');
neighber_index=find_neighber(ind_min,tri);
DT = delaunayTriangulation(xi.');
triplot(DT)
hold on
neighber_index=find_neighber(ind_min,tri);
plot(xi(1,ind_min),xi(2,ind_min),'k*','markersize',10)
plot(xi(1,neighber_index),xi(2,neighber_index),'ro','markersize',10)
for ii=1:length(neighber_index)
  delta=norm(xi(:,neighber_index(ii))-xi(:,ind_min));
%   ss(ii)=(yi(neighber_index(ii))+ymin)/2-K*delta^2/4; 
  ss(ii)=((yi(neighber_index(ii))+ymin)-2*y0)/delta^2;
end
[tt,ind_next]=min(ss);
xi=[xi (xi(:,neighber_index(ind_next))+xi(:,ind_min))/2];
yi=[yi fun((xi(:,neighber_index(ind_next))+xi(:,ind_min))/2)];
% pause
%  writeVideo(writerObj,getframe);
saveas(gcf,num2str(kk,'%0.02d'),'png')
saveas(gcf,num2str(kk,'%0.02d'),'eps')

end

%%
% DT = delaunayTriangulation(xi.');
% triplot(DT)
% hold on
% neighber_index=find_neighber(ind_min,tri);
% plot(xi(1,ind_min),xi(2,ind_min),'ks')
% plot(xi(1,neighber_index),xi(2,neighber_index),'r+','markersize',10)