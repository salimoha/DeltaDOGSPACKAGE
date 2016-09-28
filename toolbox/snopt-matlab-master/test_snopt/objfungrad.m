function [f,gradf] = objfungrad(x)
% a=100;
% b=0.0001;
% f = a*x(2)+b*sum(x.^2);
% % Gradient of the objective function:
% if nargout  > 1
%     gradf = [ 0+2*b*x(1)  ,a+2*b*x(2)];
% end

 %%
a=100;
f = a*x(2);
% Gradient of the objective function:
% % % if nargout  > 1
    gradf = [ 0,a];
% % % end
%%
% 
% a=1;
% b=1;
% f = a*x(2)+b*sum(x.^2);
% % Gradient of the objective function:
% if nargout  > 1
%     gradf = [ 0+2*b*x(1)  ,a+2*b*x(2)];
% end

%%

% 
% % 
% f = exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
% % Gradient of the objective function:
% if nargout  > 1
%     gradf = [ f + exp(x(1)) * (8*x(1) + 4*x(2)), 
%     exp(x(1))*(4*x(1)+4*x(2)+2)]';
% end