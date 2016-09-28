% test example for the ilustration and validation of nonlinear constrained
% delta dogs code
% Author: Shahrouz Alimo & Pooriya Beyhaghi 
% last update:  Oct 22nd 2015
clear all
close all
clc
global n EX KC2 
n=2;
% EX = 1;
% Method =1;

for Method =[2,1];
for EX = [1,2,3];
%%%%
rootpath = './figures/';
if ~exist(rootpath) == 1 
mkdir(rootpath);
end
% specify the objective function
fun =@(x) 1-x(2,:);
%fun=@rastriginn;
%fun=@fit8;
%fun=@plygen;
% specify the contraints functions

switch EX
    case 1
%         weierstrass
con{1}=@(x) x(2,:)-fun_weierstrass(x(1,:));
con{2}=@(x)  -x(2,:)+fun_weierstrass(x(1,:));
% bnd1 = ones(n,1)*-2; bnd2 = ones(n,1)*2;
bnd1 = [-1.7;-2]; bnd2 = [1.7;2];
case 2
%     rastriginn
con{1}=@(x) x(2,:)-rastriginn(x(1,:));
 con{2}=@(x) -x(2,:)+rastriginn(x(1,:));
case 3
 
    KC2 = 1/3;
%     KC2 = 1/3;
%     KC2 = 5;
con{1}=@(x) x(2,:)-KC2*sin(x(1,:)*pi);
con{2}=@(x) -x(2,:)+KC2*sin(x(1,:)*pi);
case 4
 
    KC2 = 3;
%     KC2 = 1/3;
%     KC2 = 5;
con{1}=@(x) x(2,:)-KC2*sin(x(1,:)*pi);
con{2}=@(x) -x(2,:)+KC2*sin(x(1,:)*pi);
end
%% constant K algorithm 
switch Method 

case 1
  Vrange = [0,0.25, 0.5, 1,1.5, 2, 5,10 ,50];
%   Vrange = [ 5,10 ,25 ,50];
%   Vrange = [1];
for Var = Vrange
    % folder
    fpath = strcat(rootpath,'/weirstrass_Method_',num2str(Method), '_const_', num2str(Var*100),'_weirstrass');
if ~exist(fpath) == 1 
mkdir(fpath);
end
    % file figure to be saved
    fName = strcat( fpath ,'/weirstrass_Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_' ,num2str(EX));
    consraints_delta_dogs_2(Var, Method, fun, con, bnd1, bnd2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%     figure_to_publish(fName, strcat('Method ', num2str(Method), ' K = ', num2str(Var)  ))
       figure_to_publish(fName)
figure(1)
    figure_to_publish(strcat(fName, '_convegance'))    
end

case 2

%% adaptive K

Vrange = [-.2,-0.1,-1e-2, 0, 1e-2, 0.1, 0.2 ];
for Var = Vrange
    % folder
    fpath = strcat(rootpath,'/weirstrass_Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
if ~exist(fpath) == 1 
mkdir(fpath);
end
if Var <0
        % file figure to be saved
    fName = strcat( fpath ,'/weirstrass_Method_',num2str(Method), '_y0_neg_', num2str(abs(Var)*100), 'ex_' ,num2str(EX));
    % file figure to be saved
else
    fName = strcat( fpath ,'/weirstrass_Method_',num2str(Method), '_y0_', num2str(Var*100), 'ex_' ,num2str(EX));
    
end
    consraints_delta_dogs_2(Var, Method, fun, con, bnd1, bnd2);
    %%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%     figure_to_publish(fName, strcat('Method ', num2str(Method), ' K = ', num2str(Var)  ))
       figure_to_publish(fName)
figure(1)
    figure_to_publish(strcat(fName, '_convegance'))    
end

end
end

end
