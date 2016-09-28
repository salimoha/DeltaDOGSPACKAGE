function [d,P,Q,sigmae2,cost,theta_P,theta_Q] = ARMAfit(x,d,theta_P,theta_Q)
% Find the best ARIMA fit for the random process x
global p q I n n1 
n=length(x);
n1=floor((n-1)/2); 
%I=abs(fft(x)).^2/n;
%I=I(2:n1+1);
t=1:n;
for i=1:n1
    landa=2*pi*i/n; 
    I(i)=abs(exp(1i*t*landa)*x)^2/(2*pi*n);
end
Lambda=2*pi*(1:n1)/n;
fun= @ whittlecost;
[A,b] = ARMA_ST_Cond(p,q);

theta=[theta_P theta_Q];
% theta= [0 -1.4595   -0.6126    1.4370    0.0418   -0.4757    0.0171    0.0541 0  -0.7232   -0.9711    0.6665    0.3230   -0.2001   -0.0372    0.0192 0] 
theta= fmincon(fun,theta,A,b);
cost=fun(theta);
theta_P=theta(1:p);
theta_Q=theta(p+1:p+q);
[P Q] = ARMA_coff(theta);
f1=dens_arima_scaled(d,P,Q);
Q1=mean(I./f1);
sigmae2=2*pi*Q1;
%figure(3)
%loglog(Lambda,sigmae2*f1/(2*pi),'-',Lambda(1:end),I(1:end),'ro','LineWidth',2)
%grid on
%xlabel('frequency','fontsize',12,'fontweight','bold')
%ylabel('epiodogram','fontsize',12,'fontweight','bold')

end

% whittle cost function
function [cost] = whittlecost(theta)
% Calculate the whittle estimate cost function for the given Hurst
% parameter
% mode 0 for FGN, mode 1 for ARIMA
global n I 
d=0;
[P Q] = ARMA_coff(theta);
f1=dens_arima_scaled(d,P,Q);
Q=mean(I./f1);
C=exp(mean(log(f1)));
cost=log(Q*C);
end
% ARIMA model coefficents
function [P Q] = ARMA_coff(theta)
%Find the polynomial representtaionation of P Q d form theta
global p q
P1=1; Q1=1;
for i=1:floor(p/2)
    P1=conv(P1,[theta(2*i) theta(2*(i-1)+1) 1]);
end
for i=1:floor(q/2)
    Q1=conv(Q1,[theta(2*i+p) theta(2*(i-1)+p+1) 1]);
end
if rem(p,2)==1
    P1=conv(P1,[theta(p) 1]);
end
if rem(q,2)==1
    Q1=conv(Q1,[theta(p+q) 1]);
end
 P1(end)=[];
 Q1(end)=[];
 P=fliplr(P1);
 Q=fliplr(Q1);

end
% The constraints for optimization
function [A,b] = ARMA_ST_Cond(p,q)
% constrcut the inequlaity ARIMA problem
A=zeros(2*(p+q),p+q); b=ones(2*(p+q),1);
for i=1:floor(p/2)
    A(4*(i-1)+1:4*i, 2*(i-1)+1:2*i)=[0 1; 0 -1; -1 -1; 1 -1];
end
for i=1:floor(q/2)
    A(2*p+4*(i-1)+1:2*p+4*i, p+2*(i-1)+1:p+2*i)=[0 1; 0 -1; -1 -1; 1 -1];
end
    if (rem(p,2)==1)
       A(2*p-1:2*p, p)=[1;-1];
    end
    if (rem(q,2)==1)
       A(2*p+2*q-1:2*p+2*q, p+q)=[1;-1];
    end

end
% calculating the scaled sepctral desity function for ARIMA model
function [f1] = dens_arima_scaled(d,P,Q)
% Calculate the scaled density function of FBM
global n n1
r=1:n1;
lamda=2*pi*r/n;
f1=(2*sin(lamda/2)).^(-2*d);
si=lamda*0+1; phi=1+lamda*0;
p=length(P); q=length(Q);
for k=1:q
    si=si+exp(1i*lamda*k)*Q(k);
end
for k=1:p
    phi=phi+exp(1i*lamda*k)*P(k);
end
f1=f1.*(abs(si).^2); f1=f1./(abs(phi).^2);


end
