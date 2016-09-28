%%

x0=[0.15; 0.15];
b=[1;1;0;0];
A = [1,0;0,1;-1,0;0,-1];
% fun=@(x)fun_constraint;
x = fmincon(@fun_constraint,x0,A,b)

%%
a=1;
x0=[0.15; 0.15];
[y0,dy0] = fun_constraint(x0);
x=x0;
for i=1:10000
    [ y1,dy] = fun_constraint(x);
      p =  -y1/dy;
      x1 = x + a*p;
      if y1-y0 <=0
      y0 = y1;
      x = x1;
%       keyboard
      else
          a = a*1/2;
      end
      
end



%%
x0=[0.15; 0.15];
[ y0,dy0] = fun_constraint(x0);
x=x0;
for i=1:100
    [ y,dy] = fun_constraint(x);
    x = x - y/dy;
end
x0
y0
x
y


%%

x0=[0.1551; 0.1552];
[ y0,dy0] = fun_constraint(x0);
x=x0;
for i=1:100
    [ y,dy] = fun_constraint(x);
    x = x - y/dy;
    
end
x0
y0
x
y


%%

