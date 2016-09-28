% test example for the ilustration and validation of nonlinear constrained
% delta dogs code
% Author: Shahrouz Alimo & Pooriya Beyhaghi 
% last update:  Oct 22nd 2015
clear all
close all
clc

rootpath = './figures/ConstCoef_1/';
if ~exist(rootpath) == 1 
mkdir(rootpath);
end
% specify the objective function
fun =@(x) 1-x(2,:);
%fun=@rastriginn;
%fun=@fit8;
%fun=@plygen;
% specify the contraints functions
global EX KC2 b
% b = pi/6;
% for Method =[1,2];
% for EX = [1,2,3];
for Method =[1,2];
for EX = [1,2,3,4];
switch EX
    case 1
         b=1/2;
    KC2 = 2;
%     KC2 = 1/3;
%     KC2 = 5;
con{1}=@(x) KC2*(x(2,:)-b*sin(x(1,:)*pi));
con{2}=@(x) KC2*(-x(2,:)+b*sin(x(1,:)*pi));
case 2
     b=1/2;
    KC2 = 5;
%     KC2 = 1/3;
%     KC2 = 5;
con{1}=@(x) KC2*(x(2,:)-b*sin(x(1,:)*pi));
con{2}=@(x) KC2*(-x(2,:)+b*sin(x(1,:)*pi));
case 3
 b=1/2;
    KC2 = 10;
%     KC2 = 1/3;
%     KC2 = 5;
con{1}=@(x) KC2*(x(2,:)-b*sin(x(1,:)*pi));
con{2}=@(x) KC2*(-x(2,:)+b*sin(x(1,:)*pi));
case 4
  b=1/2;
    KC2 = 100;
con{1}=@(x) KC2*(x(2,:)-b*sin(x(1,:)*pi));
con{2}=@(x) KC2*(-x(2,:)+b*sin(x(1,:)*pi));
end

% constant K algorithm 
switch Method 

case 1
  Vrange = [0,0.25, 0.5, 1,2, 5,10 ,50];
%   Vrange = [ 5,10 ,25 ,50];
% Vrange = [1];
for Var = Vrange
    % folder
    fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
if ~exist(fpath) == 1 
mkdir(fpath);
end
    % file figure to be saved
    fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_' ,num2str(EX));
    consraints_delta_dogs_3(Var, Method, fun, con);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%     figure_to_publish(fName, strcat('Method ', num2str(Method), ' K = ', num2str(Var)  ))
       figure_to_publish(fName)
figure(1)
    figure_to_publish(strcat(fName, '_convegance'))    
end

case 2

% adaptive K

   Vrange = [ 0,0.3,0.5,0.67,0.76 ,1];
% Vrange = [-0.5, 0, 0.3,0.5,0.67,0.76 ];
%     Vrange = [1, 0, 0.5 ];
%  Vrange = [0.5];
for Var = Vrange
    % folder
    fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
if ~exist(fpath) == 1 
mkdir(fpath);
end
if Var <0
        % file figure to be saved
    fName = strcat( fpath ,'/Method_',num2str(Method), '_y0_neg_', num2str(abs(Var)*100), 'ex_' ,num2str(EX));
    % file figure to be saved
else
    fName = strcat( fpath ,'/Method_',num2str(Method), '_y0_', num2str(Var*100), 'ex_' ,num2str(EX));
    
end
    consraints_delta_dogs_3(Var, Method, fun, con);
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
