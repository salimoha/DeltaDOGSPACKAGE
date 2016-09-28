clear all
close all
clc
L = 200;
% N = 512;
N=1024;
dx = L / N;
x = dx * (0 : N-1)';
kx = (2*pi/L) * [(0 : N/2)'; (-N/2+1 : -1)']; 
Aop = 0.5 * (kx .^ 2 - kx .^ 4);
kx(fix(N/3)+1 : N-fix(N/3)+1) = 0;
u = sin(200 * pi * x / L) + sin(340 * pi * x / L) + 0.2 * randn(N,1);
uhat = fft(u, N);
% T = 5e6;
T=1e3;
% T=600;
%  T = 5e5;

iPlot = 1;
iSave = 0;
%% CN/RKW3 Parameters

dt = 0.2;
t = 0 : dt : T;
Nt = fix(T / dt);

alphaI = [4/15 1/15 1/6];
betaI = [4/15 1/15 1/6];
betaE = [8/15 5/12 3/4];
gammaE = [0 -17/60 -5/12];
%% Simulation

% uhat_vec = zeros(N,Nt/iSave+1); cont = 1; uhat_vec(:,cont) = uhat;
% t_vec = zeros(1,Nt/iSave+1);
for k = 1 : Nt
    for rk = 1 : 3
        r = ifft(uhat, N);
        r = - r .* r / 2;
        rhat = 1i * kx .* fft(r, N);
        if rk == 1
            uhat = (uhat + dt * (betaI(rk) * Aop .* uhat + betaE(rk) * rhat)) ./ (1 - dt * alphaI(rk) * Aop);
        else
            uhat = (uhat + dt * (betaI(rk) * Aop .* uhat + betaE(rk) * rhat + gammaE(rk) * rhat_old)) ./ (1 - dt * alphaI(rk) * Aop);
        end
        if rk < 3
            rhat_old = rhat;
        end
    end
    if mod(k, iSave) == 0
%         disp(['dt = ', num2str(t(k+1))])
        cont = cont + 1;
        t_vec(cont) = t(k+1);
        uhat_vec(:,cont) = uhat;
    end
    %
%     if mod(k, iPlot) == 0
%        pause(0.1)
%        figure(1)
%        plot(x, ifft(uhat,N));
%        axis([0 L -3 3])
% %         pause(0.001)
%          drawnow
%         figure(2)
%         semilogy(kx(1 : fix(N/3)), abs(uhat(1 : fix(N/3))).^2);
%         axis([0 3 1e-8 1e-1])
% %         pause(0.001)
%  drawnow
%         figure(3)
%         loglog(kx(1 : fix(N/3)), abs(uhat(1 : fix(N/3))).^2);
%         axis([3e-2 4 1e-8 1e-1])
%         drawnow
%     end
    %
KE(k)=mean(abs(uhat(1:fix(N/3),:)).^2)/N; 
KE_1(k)=sum(abs(uhat(1:fix(N/3),:))).^2;
end
% TKE=sum(abs(uhat_vec(1:fix(N/3),:))).^2;


% save('final_long_results')
t_trans=transient_time_detector(KE);
%
figure(2);clf;
plot((5:length(KE)),KE(5:length(KE)),'k-','linewidth',2)
hold on,
grid minor
plot([(t_trans),(t_trans)],[min(KE),1],'k-.','linewidth',1.5)
xlim([5,T/dt])
ylim([0,0.9])
% figure_to_publish('../data/KUROMOTO_TRANSIENT')

%figure(1);clf;
% plot(t_vec(5:end),TKE(5:end),'k-','linewidth',2)
% hold on,
% grid minor
% plot([t_vec(t_trans),t_vec(t_trans)],[min(TKE),1],'k-.','linewidth',1.5)
% figure_to_publish('KUROMOTO_TRANSIENT')
