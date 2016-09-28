function [xi,xic] = Linear_constrained_initilazation(x0f1)
% Initialization of the linearly_cosntrained problem
% Ain x \leq bin

global n m kappa Ain bin 

% find an initial feasible point
for ii=1:m
    H{ii}=zeros(n,n);
    f{ii}=-(Ain(ii,:))';
    e{ii}=bin(ii);
end
[x0f]=feasible_point_finder(x0f1,H,f,e);
[x0r]=Initial_point_finder(x0f,H,f,e);
V = uniformsimplexrecursive(n);
bounded=1;
% prolongate it to the boundary
for ii=1:n+1
    c = (bin-Ain*x0r)./(Ain*V(:,ii)); d=Ain*V(:,ii);
    if length(d>0)==0
        printf('The feasible domain is unbounded');
        bounded==0;
        break
    else
    alpha=min(c(d>0));
   xi(:,ii)=x0r+alpha*V(:,ii);
    end
end
if bounded==0
    exit
end
err=1;
while (err>1e-3)

err=0;
for k=1:n+1
xic=xi; 
xic(:,k)=[];
xn= Lin_last_point_find(xi(:,k),xic);
err=max([err,norm(xi(:,k)-xn)]);
xi(:,k)=xn;
end
end
% find the exterior simplex
for k=1:n+1
xic(:,k)=xi*ones(n+1,1)-n*xi(:,k); 
end
O=xi*ones(n+1,1)/(n+1);
for k=1:n+1
xic(:,k)=O+kappa*(xic(:,k)-O);
end

end

