function test_FilterSQP

clc; clear all;

% x0 = [0.4,0];            % Starting guess 
% %  x0 = [0.4,1];            % Starting guess 
% % x0 = [1,1];            % Starting guess 
% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed', 'PlotFcns',{@optimplotx,...
%     @optimplotfval,@optimplotfirstorderopt} );
% % options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed', 'PlotFcns',@optimplotconstrviolation  );
% options = optimoptions(options,'GradObj','on','GradConstr','on');
% lb = [ ]; ub = [ ];   % No upper or lower bounds
% [x,fval] = fmincon(@objfungrad,x0,[],[],[],[],lb,ub,... 
%    @confungrad,options)
%%
% confungrad(x)
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','GradObj','on','GradConstr','on');
% problem.options = options;
% problem.solver = 'fmincon';
% problem.objective = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% problem.x0 = [0,0];
% problem.nonlcon = @twocone;
% x = fmincon(problem)
% 
%%
% x0=[0.4;1];
%  x0 = [0.4;-5]; % infeasible point
 x0 = [0.4;5 ]; % feasible point 
m=2;
Ain = [eye(2); -eye(2)];
% bin = [ones(4,1)];
bin = [ones(4,1)*5];
% [x,y]= spFilterSPQ(@fun, @gradfun, @hessfun, @consfun, @consgrad, x0, Ain,bin,m)
 [x,y,Xtraj,Ytraj,Ctraj]=  spFilterSPQ(@fun, @gradfun, @hessfun, @consfun, @consgrad, x0, Ain,bin,m)
%  keyboard
end
%% objective functions
function [f] = fun(x)

    f = x(2);
end

function [df] = gradfun(x)
   df = [0, 1];     
end



function d2f = hessfun(x)
    
   d2f = [0, 0;
          0, 0];
end

%% constraint functions
function [C] = consfun(x)

C{1} = -(x(2)-0.5*(x(1)-0.5).^2-0.1); % Inequality constraints
C{2} = -x(2); 
end

function [DC] = consgrad(x)
   DC{1} = [x(1)-0.5; -1];
   DC{2} = [ 0; -1]; 
%     DC= [x(1)-0.5, -1;
%          0, -1];
end
