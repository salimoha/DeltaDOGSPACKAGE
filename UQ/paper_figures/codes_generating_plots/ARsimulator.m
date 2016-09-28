function [x] = ARsimulator(coeff, simTime, x0,sigma02)
% both x0 an dp are coloun wise
p = length(coeff);
% keyboard
x=x0; coeff=fliplr(coeff);
for ii=1:simTime-p+1
   x(ii+p) = -coeff*x(ii:ii+p-1)+sqrt(sigma02)*randn; 
end

end
