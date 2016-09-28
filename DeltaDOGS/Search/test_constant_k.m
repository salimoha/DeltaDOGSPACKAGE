% test problem for constant-K-search algorithm 
clear all
clc


global n Ain bin ms m
n=2;
xc=[0.5;0.5];
n=2; ms=2; m=4;

Ain=[eye(n);-eye(n)]; bin=[ones(n,1) ;zeros(n,1)];
R2=0.5;
xi=[0 0 1 1; 0 1 0 1]; v1=[8 ;0 ;8]; v2=[0;0;1]; 
x=xc;

inter_par{1}=1; inter_par{2}= zeros(4,1); inter_par{3}= v1; inter_par{4}=xi;
inter_Par{1}{1}=1; inter_Par{1}{2}= zeros(4,1); inter_Par{1}{3}= v2; inter_Par{1}{4}=xi;
inter_Par{2}{1}=1; inter_Par{2}{2}= zeros(4,1); inter_Par{2}{3}= -v2; inter_Par{2}{4}=xi;

Search.constant=0.3651; Search.method=1;
%%
%[x y CS]=Constant_K_Search_interior(x,inter_par, inter_Par,xc,R2,Search)
%%
[x y CS]=Constant_K_Search(x,inter_par,inter_Par,xc,R2,Search);