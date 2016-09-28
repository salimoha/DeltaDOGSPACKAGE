function [H_par, f_par, d_par, Par] = QuadraticFitting(xi, yi, ymax, Res)
%UNTITLED Summary of this function goes here
%  p(x) = 1/2 x'Hx + f'x + d
% input:
% xi: the design paramter traiing set, 
% yi: function rvaluation training set,
% out put is the interpoaltion coeficients
% H: hessian of the quadratic model (H is diagonal for simplicity)
% f: linear coeficient of the model
% d: constant coeficient of the model
% N: number of data points
% n: dimension 

% weighted least square
if nargin <4
% Res = 1e2;
Res = 1e-2;
end
%N = length(yi);
% The equality constraints
[ymin, idmin] = min(yi);
xie=xi(:,yi-ymin < Res); yie=yi(yi-ymin<Res);
Aeq=[xie.^2;xie;ones(1,size(xie,2))]; 
Aeq=Aeq.'; beq=yie.';
xi(:,yi-ymin < Res)=[]; yi(yi-ymin < Res)=[];

% initialize the quadprog
[n,N]= size(xi); % xi each columns are one point
H = zeros(2*n+1,2*n+1); f=zeros(2*n+1,1);

for ii = 1:N
%% this is a modeling part for hessian so one can change 'ai' to get other models
    % for evry point in training set
    ai(:,ii) = [xi(:,ii).^2 ; xi(:,ii) ; 1];
    bi(:,ii) = -yi(ii);
    %% w4eighting strategy it can be changed to another weithting coeficients 
    w(:,ii) = 1./(yi(ii) - ymin);
    if w(:,ii) < 1e-2
        w(:,ii)=0;
    end
    H=H+w(:,ii)*ai(:,ii)*ai(:,ii)';
    f=f+w(:,ii)*ai(:,ii)*bi(:,ii);      
end

% The inequality constraints
Aineq1=-[xi.^2;xi;ones(1,size(xi,2))].'; Bineq1=-ymin*ones(size(xi,2),1);
Aineq2=[xi.^2;xi;ones(1,size(xi,2))].'; Bineq2=ymax*ones(size(xi,2),1);
% Aineq3=[xi.^2;xi;ones(1,size(xi,2))].'; Bineq3=(2*yi'-ymin).*ones(size(xi,2),1);
Aineq3=[]; Bineq3=[];
Aineq = [Aineq1; Aineq2; Aineq3];
Bineq = [Bineq1; Bineq2; Bineq3];
% Aineq*quad - Bineq
%  keyboard
% optimization scheme for quadratic programming 
% options = optimset('algorithm', 'trust-region-reflective', 'Display', 'iter');
% options = optimset('algorithm', 'interior-point-convex', 'Display', 'iter');
options = optimset('algorithm', 'active-set', 'Display', 'iter');
[quad_parameter, yy ]=quadprog(H,f,Aineq ,Bineq,Aeq,beq,[],[],zeros(2*n+1,1),options);
H_par=diag(quad_parameter(1:n));
f_par=quad_parameter(n+1:2*n);
d_par=quad_parameter(end);
%
Par.H = H;
Par.f = f;
% R100 = load('QuadCoefNew.mat')
end

function q = funquad(a0,H,f)
% for checking the algoirthm 
q =0.5*a0'*H*a0+ f'*a0;
end

