function [xi,acon] = initial_pts()
global m n bnd1 bnd2 
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points

acon = cell(m,1);
for ii=1:size(xi,2)
    ind=1:2*n;
    ind=ind(Ain*xi(:,ii)-bin>-0.01);
    for jj=1:length(ind)
    acon{ind(jj)}=[acon{ind(jj)} ii];
end
 %
% save init_pts.mat xi 
end
%
xi(1,1)=xi(1,1)+0.00001; 
end