%% paper figures
% function plot_figures()
% figure 2 example 2
close all
  figure(1); clf;
%   box on
  hold on
  %%   obj plots
    subplot(2,1,1)
%     plot(1:length(yi(1:end)),yi,'k-', 'linewidth', 3)
      plot(1:length(yi(1:end)),yi,'k-', 'linewidth', 3)
%     title('objective function')
if Ex ==1
%     y0 = 
%      plot(tt, b*sin(pi*tt), 'k--', 'linewidth', 3)
elseif Ex ==2
    x_star=[0.7,0.1];
    y0 =0.1;
    hold on
plot(1:length(yi(1:end)),  ones(length(yi(1:end)),1)*y0, 'k:', 'linewidth', 2)  
xlim([0,30])
elseif Ex==3
  x_star=[0.1559 ;0.1525]; % example island
   y_star = fun(x_star);
     hold on
     ylim([0,2])
     plot(1:length(yi(1:end)),  ones(length(yi(1:end)),1)*y_star, 'k:', 'linewidth', 2)  
 xlim([0,45])
end
    grid on
      grid minor
     set(gca,'fontsize',20)
     %% constraint plots
    subplot(2,1,2)
    plot(1:length(yi(1:end)),C{1},'k-', 'linewidth', 3)
%     title('costrained function')
if Ex ==1
%     y0 = 
%      plot(tt, b*sin(pi*tt), 'k--', 'linewidth', 3)
elseif Ex ==2
    hold on
    x_star=[0.7,0.1];
      plot(1:length(yi(1:end)), zeros(length(yi(1:end)),1), 'k:', 'linewidth', 2)  
xlim([0,30])
elseif Ex==3
  x_star=[0.1559 ;0.1525]; % example island
   y_star = fun(x_star);
     hold on
     plot(1:length(yi(1:end)),  zeros(length(yi(1:end)),1), 'k:', 'linewidth', 2)  
     xlim([0,45])
end
    grid on
      grid minor
      set(gca,'fontsize',20)
 %% points trajectory 
   figure(2); clf;
   tri=delaunayn(xi.');
%     triplot(tri,xi(1,:),xi(2,:))
    hold on
%     tt = -2:0.01:2;
if Ex ==1
     plot(tt, b*sin(pi*tt), 'k', 'linewidth', 3)
elseif Ex ==2
    x_star=[0.7,0.1];
      plot(tt, rastriginn(tt), 'k', 'linewidth', 3)  
%             plot(x_star(1),x_star(2),'kx', 'MarkerSize', 20)
            plot(x_star(1),x_star(2),'p','MarkerFaceColor','b', 'MarkerSize', 10)
    grid on
    grid minor
elseif Ex==3
       tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
x_star=[0.1559 ;0.1525]; % example island
contourf(tt,tt,U,0:0.5:0.5)
colormap bone
  grid on
  grid minor
end
   plot(x_star(1),x_star(2),'kp','MarkerFaceColor','b', 'MarkerSize', 15)   
  plot(xi(1,1:4),xi(2,1:4),'ks',  'MarkerFaceColor','w','MarkerSize', 10) 
  plot(xi(1,5:end),xi(2,5:end),'ks','MarkerFaceColor','w', 'MarkerSize', 10)

    axis square   
  
    xlim([0 1])
    ylim([0 1])
    set(gca,'fontsize',20)
    box on
%      xlim([-2 2])
%     ylim([-2 2])
    drawnow 
%    saveFigK(Search)
    
% end
%%
% for ii=1:k
% strmax = [ num2str(ii)];
% % text(xi(1,ii+4),xi(2,ii+4),strmax ,'HorizontalAlignment','right');
% text(xi(1,ii+4),xi(2,ii+4),strmax );
% end



%%
% function saveFigK(Search, Path)
% if nargin <2
if Ex ==2
    Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\example2\new_plots';
%  Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\12_week\paper_figures\example2';
elseif Ex ==3
% Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\12_week\paper_figures\island';    
% Path ='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\island';
Path ='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\island\';
end
%  Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\12_week\paper_figures\example2';
% Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\12_week\paper_figures\island';
% Path='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\example2\new_plots';
% Path ='C:\Users\shahrouz\Dropbox_aol\Dropbox\Global1_SP_1\delta_dogs\week_12\paper_figures\island';
% end
saveKon =1;
if saveKon ==1
        if Search.method ==1
%     title(strcat('Position of the points',' for K = ', num2str(Search.constant), '  KC2 = ', num2str(KC2) ))
    figure(2)   
    saveas(gcf,strcat(Path,'/K_', num2str(Search.constant), 'iter_',num2str(k) ), 'png')
    saveas(gcf,strcat(Path,'/K_', num2str(Search.constant),'iter_',num2str(k) ), 'epsc2')
    saveas(gcf,strcat(Path,'/K_', num2str(Search.constant), 'iter_',num2str(k) ), 'fig')
    figure(1)   
    saveas(gcf,strcat(Path,'/fun_con_K_', num2str(Search.constant), 'iter_',num2str(k) ), 'png')
    saveas(gcf,strcat(Path,'/fun_con_K_', num2str(Search.constant),'iter_',num2str(k) ), 'epsc2')
    saveas(gcf,strcat(Path,'/fun_con_K_', num2str(Search.constant), 'iter_',num2str(k) ), 'fig')
        elseif Search.method ==2
%    title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant), '  KC2 = ', num2str(KC2) ))
    figure(2)
    saveas(gcf,strcat(Path,'/y0_', num2str(ceil(Search.constant*100)), 'iter_',num2str(k) ), 'png')
    saveas(gcf,strcat(Path,'/y0_', num2str(ceil(Search.constant*100)),'iter_',num2str(k) ), 'epsc2')
    saveas(gcf,strcat(Path,'/y0_', num2str(ceil(Search.constant*100)), 'iter_',num2str(k) ), 'fig')
    figure(1)
    saveas(gcf,strcat(Path,'/fun_con_y0_', num2str(ceil(Search.constant*100)), 'iter_',num2str(k) ), 'png')
    saveas(gcf,strcat(Path,'/fun_con_y0_', num2str(ceil(Search.constant*100)),'iter_',num2str(k) ), 'epsc2')
    saveas(gcf,strcat(Path,'/fun_con_y0_',num2str(ceil(Search.constant*100)), 'iter_',num2str(k) ), 'fig')
        end
end
        
        
% end
