function plot_figures_3(xi,yi,C, tt, Ex)
if nargin < 3 
    tt = 0:0.1:1;
end
if nargin <5
    Ex=3;
end
global Search KC2  b

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
%     tt = -2:0.01:2;

if Ex ==1
     plot(tt, b*sin(pi*tt), 'k', 'linewidth', 3)
elseif Ex ==2
      plot(tt, rastriginn(tt), 'k', 'linewidth', 3)  
elseif Ex==3
       tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.5:0.5)
end
  
%   

    plot(xi(1,:),xi(2,:),'s')
    
    if Search.method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant), '  KC2 = ', num2str(KC2) ))
    elseif Search.method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant), '  KC2 = ', num2str(KC2) ))
    end
    xlim([0 1])
    ylim([0 1])
%      xlim([-2 2])
%     ylim([-2 2])
    drawnow  


end
