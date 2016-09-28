function [mu,sigma,transtime]= Alpha_Dogs_Functionevaluation(x,index,new,transtime,a1,a2,b1,b2,ss0,rk)
% The function evaluation process for alpha dogs.
% There are two possible scenarios: a) New point b) Existing point

% Inputs
% x: parameters of the desired point
% index: the index of the point
% new: a Boolian: new point=1 , existingpoint=0;
% trasntime

% outputs
% mu: Best estimate at this point for the cost function
% sigma: The uncertainty function at this point
% transtime: The transient time for this point

% The parameter real problem=1 is fo hydofoil, and real problem==1 is for
% test problem

mainfolder=pwd;

real_problem=1;
if real_problem==0
    x0=[0.6;0.7];
    fun=@(x) norm(x-x0)^2;
    if new==1
       transtime=0; sigma=1; 
       mu=fun(x)+randn*sigma;
    else
       transtime=0; 
       mu1=fun(x)+randn;
       mu=(mu1*sigmaT(index)^2+yiT(ind))/(sigmaT(index)^2+1);
       sigma=sqrt(1/sigmaT(index)^2+1); 
    end
else    
if new==1
   % Calculate the actual variables
x=ss0+rk.*x;
x3=(b1-a1(1:2)*x(1:2))/a1(3);
x6=(b2-a2(1:2)*x(3:4))/a2(3);
x=[x(1:2) ;x3 ;x(3:4) ;x6 ;x(5:8)];
    Geometry_calculation(x', index, mainfolder);
end  


mu=0;
sigma=0;
transtime=0;
% Calculate the text files and initial time series: ?????

% run NFA: 

% Read the time series file
%cd([mainfolder '\' num2str(1e8+index)]);
%a=import('effiencyseries'); timeseries=a.data;
% Calculate the statitical analysis
%if new==1
%Initial_mod(1)=0; Initial_mod(2)=0;
%else
%Initial_mod(1)=1; Initial_mod(2)=transtime; 
%end

%[mu,sigma,Initial_mod] = alpha_statitical_analyze(timeseries,Initial_mod);
%transtime=Initial_mod(2);
end
end


function eff = Geometry_calculation(params, index, mainfolder)

% take the paramters of the gemotery and calculate the geometry

% generate the new folder for this file
a=num2str(1e8+index); mkdir(a);
%cd([mainfolder '/Geometry_generator'])  % go to the geometry constrcution format
params=[params 0 0 0 0];   % Add the additional parameters
pyparams = py.list(params); % change the parameter for the python solver
Nchord = 12;   Nspan = 101; % the number of disceretizations in chord and span
py.plygen.preAVL(pyparams,Nchord,Nspan); % generate the geometry                                                  % preprocessing: writing the AVL input files                                                          % convert the output to a python list

%keyboard
copyfile('foil.avl',[mainfolder '/' a])
copyfile('foil.ply',[mainfolder '/' a])

copyfile('foil.polar',[mainfolder '/' a])
copyfile('avl',[mainfolder '/' a])
copyfile('avl_script',[mainfolder '/' a])
copyfile('bezier.py',[mainfolder '/' a])
copyfile('bezier.pyc',[mainfolder '/' a])
copyfile('plygen.py',[mainfolder '/' a])
copyfile('plygen.pyc',[mainfolder '/' a])


end


function [mu,sigmaT,Initial_mod] = alpha_statitical_analyze(x,Initial_mod)
% Statistical analysis of the time_series vector x.
% mu: estimate of the infinite_time average.
% sigma: uncertainty of the estimate mu.


% Delete the initial transient part
figure(1)
plot(x)
if Initial_mod(1)==0
t= transientpart(x);
Initial_mod(2)=t;
else
t=Initial_mod(2);
end

x=x(t+1:end);
% Normalize x
mu=mean(x); S=std(x); 
x=(x-mu)./S;

% estimate sigmaT of the Final point 
sigmaT=S*stationary_model(x);
% 

end

function t= transientpart(x)
N=length(x); ns=50;
for i=1:1.8*ns
    s(i)=std(x(floor(N/(2*ns))*i:end));
    n(i)=length(x(floor(N/(2*ns))*i:end));
end
figure(2)
plot(floor(N/(2*ns))*(1:1.8*ns),(s./sqrt(n))/std(s./sqrt(n)))
inter_par=interpolateparametarization((1:ns)/ns,(s./sqrt(n))/std(s./sqrt(n)),1);
[t,ind]=min(s./sqrt(n)); tx=inter_min(ind/ns, inter_par,0);
t=floor(tx*N/2);
end

function sigmaT=stationary_model(x)
N=length(x);
sigma=std(x); x=x/std(x);
p=2; q=1; d=0;

% Estimate the integral time scale
rho=autocorr(x,1000);
f= fit((1:length(rho))',rho','exp1');
tauf=-1/f.b;
T=round(sqrt(1/2*tauf)); T=min(T,floor(N/10));

% Estimate sigma(T)
for i=1:N-T-1 yy(i)=mean(x(i+1:i+T)); end
sigmaT=std(yy)*sqrt(T/N);
% put safety factor for small data
if T==floor(N/10)
    sigmaT=2*sigmaT;
end
end


