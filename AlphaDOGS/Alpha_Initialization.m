% The code to Ininialize the Alpha-Dogs
clear all
close all
clc

warning off 
global n m Ain bin ss0 rk

% Dimension
n=8; 

cd FilesToPath
addTopath;
cd ..

setenvAVL
% The upper bound and lowe bound for each variable
%ss0=[0.01 1.0472 0.01 0.5 0.1 0.7 0.3 0 0 0].';
%ss1=[pi/2 pi/2 pi/2 1.5 0.5 1.5 0.5 1 1 1].';
%rk=ss1-ss0;

% Modify the equality constraints for dimension reduction 

% Initilaize Ain, Bin
Ain=[]; bin=[]; ss0=[]; rk=[];

% First three paramters
a1=[1 1 1]; b1=pi/2; % equality constraints for them a^T x=b
lb1=[0.01 1.0472 0.01]; ub1=[pi/2 pi/2 pi/2];



[A,B,s,r]=actual_scaling_lincon([],[],ub1',lb1',a1,b1,3,1);
A=[A zeros(size(A,1),n-2)];
Ain=[Ain; A]; bin=[bin; B];
ss0=[ss0 ;s]; rk=[rk ;r];


% Second three paramteres
a2=[1 1 1]; b2=2; % equality constraints for them a^T x=b
lb2=[0.01 0.1 0.7]; ub2=[1.5 0.5 1.5];


[A,B,s,r]=actual_scaling_lincon([],[],ub2',lb2',a2,b2,3,1);
A=[zeros(size(A,1),2) A zeros(size(A,1),n-4)]; 
Ain=[Ain; A]; bin=[bin; B];
ss0=[ss0 ;s]; rk=[rk ;r];

% Last Four parameters
lb3=[0.3 0 0 0]; ub3=[0.5 1 1 1];

s=lb3'; r=(ub3-lb3)';
Ain=[Ain;[zeros(4,4) eye(n-4)];[zeros(4,4) -eye(n-4)]]; bin=[bin; ones(n-4,1); zeros(n-4,1)];
ss0=[ss0 ;s]; rk=[rk ;r];
m=size(Ain,1); 
% Calculate the position of initial points
  [xie ,xic]= Linear_constrained_initilazation(Ain,bin);
   O=xie*ones(n+1,1)/(n+1);
   xie=[O xie]; xiT=[xic xie];
   [xiT ,xie]=inter_add_fix(xiT,xie);
   
% Calculate at the initial points
for ii=n+2:size(xiT,2)
[yiT(:,ii), sigmaT(ii),transtime(ii)]=Alpha_Dogs_Functionevaluation(xiT(:,ii),ii,1,0,a1,a2,b1,b2,ss0,rk);
end
yiT(1:n+1)=inf; sigmaT(1:n+1)=inf; transtime(1:n+1)=inf;
