% test example for the ilustration and validation of nonlinear constrained
% delta dogs code
% Author: Shahrouz Alimo & Pooriya Beyhaghi 
% last update:  Oct 22nd 2015
clear all
close all
clc
%%
% Clear variables; be more agressive if the user has specified a preference to do so
clear cleanupTask; clear % important to clear the cleanupTask variable before everything else
if( ispref('si','clear') && getpref('si','clear') ), clear all; close all; clc; end

%%
% [ conf ] = readConf( conf_path , convertFlag )



%%
rootpath = './figures/Example_RAS1/';
if ~exist(rootpath) == 1 
mkdir(rootpath);
end
% specify the objective function
% fun =@(x) 1-x(2,:);
%fun=@rastriginn;
%fun=@fit8;
%fun=@plygen;
% specify the contraints functions
global EX KC2 b
% b = pi/6;
% for Method =[1,2];
% for EX = [1,2,3];
% EX = 'RASTRIGINN_1'

% for Method =[1,2];
% for EX = [1,2,3,4];
for Method =[1,2];
for EX = [1:4];
switch EX
    case 1
    b=1/2;
    KC2 = 1;
% specify the objective function
fun =@(x)x(2,:);
% specify the contraints functions
con{1}=@(x) KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
case 2
     b=1/2;
    KC2 = 5;
% specify the objective function
fun =@(x)x(2,:);
% specify the contraints functions
con{1}=@(x) KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
case 3
 b=1/2;
    KC2 = 10;
% specify the objective function
fun =@(x)x(2,:);
% specify the contraints functions
con{1}=@(x) KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
case 4
  b=1/2;
    KC2 = 100;
% specify the objective function
fun =@(x)x(2,:);
% specify the contraints functions
con{1}=@(x) KC2*(x(2,:)-rastriginn(x(1,:)));
con{2}= @(x) -con{1}(x);
end

% constant K algorithm 
switch Method 
case 1
  Vrange = [0,0.25, 0.5, 1,2, 5,10 ];
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
    consraints_delta_dogs(Var, Method, fun, con);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%     figure_to_publish(fName, strcat('Method ', num2str(Method), ' K = ', num2str(Var)  ))
       figure_to_publish(fName)
figure(1)
    figure_to_publish(strcat(fName, '_convegance'))    
end

case 2

% adaptive K

%    Vrange = [ 0,0.3,0.5,0.67,0.76 ,1];
Vrange = [-0.5, 0, 0.3,0.5,0.67,0.76 ];
%     Vrange = [1, 0, 0.5 ];
%  Vrange = [0.5];
%  Vrange = [0.];
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
    consraints_delta_dogs(Var, Method, fun, con);
    %%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
       figure_to_publish(fName)
figure(3)
       figure_to_publish(strcat(fName, '_c_i_interpolation_at_0'))
           
figure(1)
    figure_to_publish(strcat(fName, '_convegance'))    
end

end

end
end
