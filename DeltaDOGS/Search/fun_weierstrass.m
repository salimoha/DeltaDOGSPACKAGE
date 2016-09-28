function y = fun_weierstrass(x, N)
% weierstrass function 1D example in delta dogs paper
x = -pi/2+(0.7+pi/2)*x;
if nargin <2
    N = 300;
end
y=0;
for i=0:N
    y = y + 1/2^i * cos(3^i*pi*x);
end
y=1/4*(y+2); 
end