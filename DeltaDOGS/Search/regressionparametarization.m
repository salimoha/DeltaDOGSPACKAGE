function [inter_par,yp]= regressionparametarization(xi,yi,sigma,inter_method)
global n beta

rho0=10;
rho=rho0; c=0.8;

if inter_method==1
 while 1
    N = size(xi,2); A = zeros(N,N);
  % calculate regular A matrix for polyharmonic spline
for ii = 1 : 1 : N
    for jj = 1 : 1 : N
        A(ii,jj) = ((xi(:,ii) - xi(:,jj))' * (xi(:,ii) - xi(:,jj))) ^ (3 / 2);
    end
end
% Add the noise term
for ii=1:N
    A(ii,ii)=A(ii,ii)+rho*sigma(ii)^2;
end

V = [ones(1,N); xi];
A = [A V'; V zeros(n+1,n+1)];
wv = A \ [yi.'; zeros(n+1,1)]; % solve the associated linear system
inter_par{1}=1;
inter_par{2} = wv(1:N); inter_par{3} = wv(N+1:N+n+1); 
inter_par{4}= xi;

for ii=1:N
    yp(ii)=interpolate_val(xi(:,ii),inter_par);
end

if max(abs(yp-yi)./sigma)<beta
   break
else
    rho=c*rho;
end
 end
if inter_method==2
    inter_par=[]; yp=[];
 disp(sprintf( 'Wrong Interpolation method') );
end
end