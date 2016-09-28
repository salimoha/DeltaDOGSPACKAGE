close all; clc;
% xi=xE; yi=yE;
% load 
fs=14; ls=1.5;
figure(1);clf;
% C0=C; yi0=yi;yi=yi-FUN.upb;
% [  C ] = const_violation( C, ConBound );
for i = 1:length(yi)
%      [VM0(i),indm(i)] = max([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i),C{7}(i) ]);
    [VM0(i),indm(i)] = max([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i) ]);
    % VM(i) = max([yi(i),C{1}(i), C{3}(i) ]);
%     [Vm0(i),ind(i)] = min([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i),C{7}(i) ]);
%     [Vm0(i),ind(i)] = min([yi(i),C{1}(i),C{2}(i), C{3}(i),C{4}(i),C{5}(i), C{6}(i) ]);
%     VM0(i) = max([yi(i),C{1}(i), C{3}(i) ]);
end
% VM0 = VM;
subplot(3,1,1); 
hold on
plot(VM0,'k','linewidth',ls)
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',1.5)
% plot(Vm0,'k--','linewidth',3)
% plot(VM,'b--','linewidth',3)
hold off
grid on
grid minor
 ylim([-0.1,2])
set(gca,   'FontSize', fs);



subplot(3,1,3); 
hold on
% trajectory of maximum violation
for ii = 1:length(yi)
    [ym(ii),indx(ii)] = min(VM0(1:ii));
    Xvm(:,ii) = xi(:,indx(ii));
end
hold on
plot(Xvm(1,:),'--r','linewidth',1.5)
plot(Xvm(2,:),'-k','linewidth',1.5)
plot(Xvm(3,:),':b','linewidth',1.5)

% plot(VM,'k--','linewidth',3)
hold off
grid on
grid minor
ylim([-0.1,1.1])
% title('best point trajectory')
% set(ah, 'FontSize', op.fontsize);
set(gca,   'FontSize', fs);




subplot(3,1,2); 
hold on
% trajectory of maximum violation
plot(xi(1,:),'--r','linewidth',1.5)
plot(xi(2,:),'-k','linewidth',1.5)
plot(xi(3,:),':b','linewidth',1.5)
% plot(xi(1,:),'k-','linewidth',3)
% plot(xi(2,:),'b-.','linewidth',3)
% plot(xi(3,:),'r--','linewidth',3)
% plot(VM,'k--','linewidth',3)
hold off
grid on 
grid minor
ylim([-0.1,1.1])
set(gca,   'FontSize', fs);
% samexaxis('abc','xmt','on','ytac','join','yld',1)
   samexaxis('XLim',[0,length(C{1})+10],'YAxisLocation','alternate2' , ...
      'abc','xmt','on','ytac','join','yld',1., 'YTickAntiClash',1)
% cmfile = './ci_new';
cmfile = './ci_basic';
 figure_to_publish(cmfile)

%              saveas(gcf,cmfile, 'fig');       
% %             saveas(gcf,cmfile, 'png');
%              saveas(gcf,cmfile, 'epsc2'); 
% 
% 
% 
%   subplot(3,1,1); plot(randn(100,1),randn(100,1),'x'); ylabel('QF')
%   subplot(3,1,2); plot([-1 0 .5 1],[0 100 100 10],'x-'); ylabel('HT');
%   subplot(3,1,3); plot(randn(100,1),randn(100,1)*33,'x'); ylabel('DV');
%%

figure(2);clf;
subplot(3,1,1);
plot(C{1},'k-','linewidth',ls); 
% ylabel('\sigma_{\infty}', 'Fontsize',fs)
hold on
% plot(C{2},'--','linewidth',3)
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',ls)
 grid on
 grid minor
 set(gca,   'FontSize', fs);
ylim([-2,2])
 
 subplot(3,1,2);  
hold on
 plot(yi,'k','linewidth',ls);
%  ylabel('A^{(4)}', 'Fontsize',fs)
ylim([-.1,2])
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',ls)
 grid on
 grid minor
 set(gca,   'FontSize', fs);
 subplot(3,1,3);
 hold on
 %
 
plot(C{3},'k-','linewidth',ls); 
% ylabel('\delta', 'Fontsize',fs)
plot([0,length(C{1})] , [0,0], 'r-.','linewidth',ls)
hold off
grid on
grid minor
ylim([-2,2])
set(gca,   'FontSize', fs);
% samexaxis('XLim',[0,90],'XMinorTick','abc','xmt','on','ytac','join','yld',1)
%  samexaxis('XLim',[0,90],'YAxisLocation','alternate2' , 'YLabelDistance', 5, ...
%      'abc','xmt','on','ytac','join','yld',1., 'YTickAntiClash',1)
  samexaxis('XLim',[0,length(C{1})+10],'YAxisLocation','alternate2' , ...
     'abc','xmt','on','ytac','join','yld',1., 'YTickAntiClash',1)

% cmfile = './cv3_new';
cmfile = './cv3_basic';
figure_to_publish(cmfile)

%              saveas(gcf,cmfile, 'fig');       
% %             saveas(gcf,cmfile, 'png');
%              saveas(gcf,cmfile, 'epsc2'); 