% In this code we test the nonlinear delta-Dogs algorithm
%
% 
%
%
%%
clear all
close all
clc
global  bnd1 bnd2 kappa rk ss0 y00 
global n K y0 mc ms fun Ain bin
kappa=1.2;
% setenvAVL
alpha=1;
glo=0;
y0 = -1
%% specify the objective function
%fun=@alisonfit;
% fun=@fit8;
%fun=@plygen;
%fun=@rastriginn;

% parameters
n=8; mc=0; meq=0; Np=100; Nr=2; pho=1; y_max=3; y00=0.1; ms=0;
% K=1000;
% K=0;
K=2;
% interpolaion strategy
inter_method=1;
% upper and lower bouds
ss0=[10 1e-4 10 1e-4 10 1e-4 100 1e-4]';
ss1=[1000 1 100 2.5 20 1.6 10000 pi/2]';

% the inequality constraints
Ain=[]; bin=[];
% the equlaity constraints
Aeq=[]; beq=[];
% initial calculations
bnd1=zeros(n,1); bnd2=ones(n,1); 
rk=ss1-ss0;
%% input the inequality constraints
%Ain=[]; bin=[];
% if mc~=0
% bin=bin-Ain*ss0;
% Ain=Ain.*repmat(rk.',mc,1);
% for ii=1:mc
%     l=norm(Ain(ii,:));
%      Ain(ii,:)=Ain(ii,:)/l;
%      bin(ii)=bin(ii)/l;
% end
% end
% mc=2*n+mc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% test example
Xi = [0, 0, 3; 0, 3, 0];
% x=[0, 0, 2]';
% y=[0, 2, 0]';
fun =@(x) sum((1-x).^2);
% state cinstrant
c1 =@(x) x(1,:).^2 - x(2,:);
c2 = @(x) -(x(1,:).^2 - x(2,:));
ms=2;
Yc1 = c1(Xi)
Yc2 = c2(Xi)
% Xi * ones(N,1)/(N+1)
X0 = mean(Xi,2);
% X0=[0.5;0.5];
% X0=[1.5;0.5];
y = fun(Xi);
Ain = [ 1,0; 0,1;-1,0;0,-1]; 
bin = [2;2;0;0];
%%
inter_par = interpolateparametarization(Xi,y,1)
inter_Par{1}= interpolateparametarization(Xi,Yc1,1)
inter_Par{2}= interpolateparametarization(Xi,Yc2,1)
%%
Par.inter_par_S = inter_par;
Par.inter_par_C = inter_Par;
N = size(Xi,1);
[xc, R2] = circhyp(Xi, N);
%% input
x = X0;
[x y p  cse]=Constant_K_Search(X0,Par,xc,R2)

%% checking KKT condition
% 1st order condition
ys = costSearch(x,Par.inter_par_S, xc,R2)
gs = kgradSearch(x,Par.inter_par_S, xc,R2)
H = khessSearch(x,Par.inter_par_S, xc,R2)
left = H*p 
right = -gs

%

CSl =[]; Jl=[]; %Hc=[];
for l=1:ms
CSl=vertcat(CSl,costSearch(x,inter_Par{l},xc, R2));
Jl=vertcat(Jl, kgradSearch(x,inter_Par{l},xc, R2)');
end

left_c = Jl*x
right_c  = -CSl
A = vertcat(Ain, Jl);
b = vertcat(bin - Ain*x, -1.*CSl);

% 
% %% from the function Constant K Search
% % global n Ain bin m ms mc
% % parameters of backtracking
% cc=0.01;
% rho=0.9; 
% % Initialize the point in the simplex
% iter=1; 
% % Calculate the Newton direction   
% % search function 
% y=costSearch(x,inter_par,xc, R2)
% % y=cost(x,inter_par,xc, R2);
% g=kgradSearch(x,inter_par,xc,R2)
% H=khessSearch(x,inter_par,xc,R2)
% H=modichol(H,0.1,20);
% H=(H+H')/2;
% % constrained function
% CSl =[]; Jl=[]; %Hc=[];
% for l=1:ms
% CSl=vertcat(CSl,costSearch(x,inter_Par{l},xc, R2));
% Jl=vertcat(Jl, kgradSearch(x,inter_Par{l},xc, R2)');
% % Hcl=khessSearch(x,inter_Par{l},xc, R2);
% % Hcl=modichol(Hcl,0.1,20);
% % Hcl=(Hcl+Hcl')/2;
% % Hc=vertcat(Hc, Hcl);
% end
% %% linearized constrained for sqp
% A = vertcat(Ain, Jl);
% b = vertcat(bin - Ain*x, -1.*CSl);
% %%
% options=optimoptions('quadprog','Display','none');
% p=quadprog(H,g,A,b,[],[],[],[],zeros(n,1),options);
% % Perform the hessian modification
% %%
% 
% %%
% 
% %%
