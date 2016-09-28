% Alpha_Dogs for Linearly constrained problems
% Specifically designed for Hydrofoil Problem
% Author: Pooriya Beyhaghi


% AddToPath and seten
cd FilesToPath
addTopath;
cd ..


global n beta m Ain bin   
% n: Dimension, min:number of additional inequality constraints 

%n=10; min=0;  
n=2; min=0;

% Number_data_max:Maximum number of data points that has to be used.
% Number_data_max:Maximum computational Time at a point has to used.

Number_data_max=100; T_max=100;

% The upper and lower bound for variables
%ss0=[10 1e-4 10 1e-4 10 1e-4 100 1e-4]';
%ss1=[1000 1 100 2.5 20 1.6 10000 pi/2]';
ss0=[0;0]; ss1=[1;1];


% The inequality constraints Ain x \leq bin
Ain=[]; bin=[];

% Regrator parameter
beta=2;

% Scaling and changing the variables

% Initial calculations
bnd1=zeros(n,1); bnd2=ones(n,1); 
rk=ss1-ss0;

if min~=0
% normalize the additional inequlaity constraints
bin=bin-Ain*ss0;
Ain=Ain.*repmat(rk.',m,1);
for ii=1:m
    l=norm(Ain(ii,:));
     Ain(ii,:)=Ain(ii,:)/l;
     bin(ii)=bin(ii)/l;
end
end

m=2*n+min;

% Input the equality constraints
% if meq~=0
% beq=beq-Aeq*ss0;
% Aeq=Aeq.*repmat(rk.',meq,1);
% x0=Aeq\beq;
% V=null(Aeq);
% else
%     x0=zeros(n,1); V=eye(n);
% end


Ain=[Ain; eye(n);-eye(n)];
bin=[bin; bnd2 ;-bnd1];
%bin=bin-Ain*x0; Ain=Ain*V;
%n=n-meq;

% Calculate the position of initial points
  [xie ,xic]= Linear_constrained_initilazation(0.1*ones(n,1));
   O=xie*ones(n+1,1)/(n+1);
   xie=[O xie]; xiT=[xic xie];
   [xiT ,xie]=inter_add_fix(xiT,xie);
   
% Initial calculation ??????
yiT=(1:size(xiT,2))*0; sigmaT=yiT; transtime=yiT;
for ii=n+2:size(xiT,2)
    %xmr=ss0+rk.*(x0+V*xiT(:,ii));
    xmr=ss0+rk.*xiT(:,ii);
[yiT(:,ii), sigmaT(ii),transtime(ii)]=Alpha_Dogs_Functionevaluation(xmr,ii,1,0,xiT,yiT,sigmaT);
end
yiT(1:n+1)=Inf; sigmaT(1:n+1)=Inf; transtime(1:n+1)=Inf;

% Perform the iterations
iter_max=0; inter_method=1; y0=0;
for k=1:1

% Find the regression
   [inter_par,yp]=regressionparametarization(xiT(:,n+2:end),yiT(n+2:end),sigmaT(:,n+2:end),inter_method);
   
% Disceret search function
   sd=(min(yp,2*yiT(n+2:end)-yp)-y0)./sigmaT(n+2:end);
   [t,ind_min]=min(yiT); [t,ind_s]=min(sd);
   if ind_s~=ind_min
        % probably changed: Improve the measurement at xiT(:,ind_s)
        [yiT(:,ind_s), sigmaT(ind_s)]=improve(ind_s, transtime(ind_s));
        % Improvement 
   else
       [y,ind]=min(yiT); xmin=xiT(:,ind); [xm ym]=inter_min(xmin,inter_par); 
    if (ym>y0) 
        [t,ind]=min(yiT);
        [xm cse]= tringulation_search_bound(inter_par,xiT,tri);
    end
      deltam=mindis(xm,xiT);
   if sigmaT(ind_min)<deltam*dxdt
    % scaled point
       xmr=ss0+rk.*(x0+V*xm);
       [yim, sigmam, transm]=initial_cal(xmr,size(xiT,2)+1); 
   %  Modify the data set
   sigmaT=[sigmaT sigmam]; yiT=[yiT yim]; transtime=[transtime transm];
   else
       xmr=ss0+rk.*(x0+V*xiT(:,ind_min));
       [yiT(:,ind_min), sigmaT(ind_min)]=improve(ind_min,transtime(ind_min));    
   end
   end
   
   % plotting 
   figure(1)
    subplot(2,1,1)
    plot(1:length(yiT(n+2:end)),yiT(n+2:end),sigmaT(n+2:end),'.')
    ylim([0 2])
    subplot(2,1,2)
    plot(xiT(:,n+2:end)')
    drawnow
    
  
end