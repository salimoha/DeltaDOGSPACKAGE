function [c,ceq,DC,DCeq] = confungrad(x)
c(1,1) = -(x(2)-0.5*(x(1)-0.5).^2-0.1); % Inequality constraints
c(2,1) = -x(2); 
% No nonlinear equality constraints
ceq=[];
% Gradient of the constraints:
% % % if nargout > 2
    DC= [x(1)-0.5, -1;
         0, -1];
    DCeq = [];
    % for snopt
%     c=-c;
%     DC = -DC;
% % % end

%%
% c(1,1) = 1.5 + x(1) * x(2) - x(1) - x(2); % Inequality constraints
% c(2,1) = -x(1) * x(2)-10; 
% % No nonlinear equality constraints
% ceq=[];
% % Gradient of the constraints:
% if nargout > 2
%     DC= [x(2)-1, -x(2);
%         x(1)-1, -x(1)];
%     DCeq = [];
% end