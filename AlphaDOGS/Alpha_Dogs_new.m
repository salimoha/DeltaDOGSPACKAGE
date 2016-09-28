%for n=2:5
global n lb ub Ain bin

n=2;
sigma0=0.3;
%fun=@(x) 5*norm((x-0.3)).^2+sigma0*randn;
%funr=@(x) 5*norm((x-0.3)).^2;
fun=@(x) -sum(500*x.*sin(sqrt(abs(500*x))))/250+sigma0*randn;
funr=@(x) -sum(500*x.*sin(sqrt(abs(500*x))))/250;
lb=0*ones(n,1); ub=ones(n,1);
Ain=[eye(n) ;-eye(n)]; bin=[ones(n,1); zeros(n,1)];
for ff=1:3
    clearvars -except ff n lb ub fun regret estimate mesh datalength sigma0 funr sigma0
    clc

plot_index=1;
% maximum number of iterations:
iter_max=300;
%iter_max=10;

% interpolation strategy:
inter_method=1;

% Calculate the Initial trinagulation points
Nm=8; L0=1; K=2; 
xE=rand(n,n+2);
xE=round(xE*Nm)/Nm;

% Calculate the function at initial points
for ii=1:size(xE,2)
    yE(ii)=fun(lb+(ub-lb).*xE(:,ii));
    T(ii)=1;
end
xU= bounds(zeros(n,1),ones(n,1),n);
%K0=max(K0,0.1);
% initialize Nm, L, K
L=L0;
for k=1:iter_max
 %   keyboard
   [inter_par,yp]=regressionparametarization(xE,yE,sigma0./sqrt(T),inter_method);
   %inter_par=interpolateparametarization(xi,yi,inter_method);
  K0=range(yE);
   % Calculate the discrete function. 
   [yt,ind_out]=min(yp+sigma0./sqrt(T));
   sd=min(yp,2*yE-yp)-L*sigma0./sqrt(T);
   [ypmin,ind_min]=min(yp);
   [yd,ind_exist]=min(sd); xd=xE(:,ind_exist);
   %ind_min=1; ind_exist=1;
   %yd=inf;
   if (ind_min~= ind_min )%|| ind_min~= ind_out)
    yE(ind_exist)=((fun(xd))+yE(ind_exist)*T(ind_exist))/(T(ind_exist)+1);
    T(ind_exist)=T(ind_exist)+1;
   else
      
     % if  sigma0./sqrt(T(ind_exist))<0.01*range(yE)*(max(ub-lb))/Nm
     %     yd=inf;
     % end
    % Calcuate the unevaluated function:
    clear yu
    yu=[];
    if size(xU,2)~=0;
    for ii=1:size(xU,2)
               yu(ii)=(interpolate_val(xU(:,ii),inter_par)-min(yp))/mindis(xU(:,ii),xE);
    end
    end
   if (size(xU,2)~=0 && min(yu)<0)
       [t,ind]=min(yu);
       xc=xU(:,ind); yc=-inf;
       xU(:,ind)=[];
   else
   % Minimize s_c^k(x)
       while 1
            [xc yc] = tringulation_search_bound_constantK(inter_par,[xE xU],K*K0,ind_min);
           if interpolate_val(xc,inter_par)<min(yp)
               xc=round(xc*Nm)/Nm;
                break
           else
               xc=round(xc.*Nm)/Nm;
               if mindis(xc,xE)<1e-6
                   break
               end
               [xc,xE,xU,success]=points_neighbers_find(xc,xE,xU); 
               if success==1
                      break
               else
                      yu=[yu (interpolate_val(xc,inter_par)-min(yp))/mindis(xc,xE)]; 
               end
           end
       end
  if size(xU,2)~=0 
     if min(yu)<(interpolate_val(xc,inter_par)-min(yp))/mindis(xc,xE)
         [t,ind]=min(yu);
         xc=xU(:,ind); yc=-inf; 
         xU(:,ind)=[];
     end
   end
   end
   if (Nm>128 && sigma0./sqrt(max(T))<range(yE)*(max(ub-lb))/Nm)
       break
   end
    
    % Minimize S_d^k(x)
    if (mindis(xc,xE)<1e-6)
        %if (max(T)>200)
        K=2*K; Nm=2*Nm; L=L+L0;
        %else
        %    yd=-inf;
        %end
    end
        
    if yc<yd
        %keyboard
        if mindis(xc,xE)>1e-6
            xE=[xE xc]; yE=[yE fun(lb+(ub-lb).*xc)]; T=[T 1];
        end
    else
         yE(ind_exist)=((fun(lb+(ub-lb).*xd))+yE(ind_exist)*T(ind_exist))/(T(ind_exist)+1);
         T(ind_exist)=T(ind_exist)+1;
    end
   end
  regret(ff,k)=funr(xE(:,ind_min));
  estimate(ff,k)=yE(ind_min);
  datalength(ff,k)=length(xE);
  mesh(ff,k)=Nm;
  %if Nm>60
  %    break
  %end

  if plot_index
    subplot(2,1,1)
    errorbar(1:length(yE),yE,sigma0./sqrt(T))
    subplot(2,1,2)
    plot(xE.')
    drawnow
  end
end    
end

if 0
save('alpha_shwefel_n==3')
figure(1)
plot(estimate.','--','linewidth',1)
figure(2)
plot(regret.','-','linewidth',1)
figure(3)
plot(datalength.','-','linewidth',1)
end
%save(['Alpha_parab' num2str(n)])
%end
  

