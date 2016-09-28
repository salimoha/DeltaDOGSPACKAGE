function [x,fval,INFO]=test_snopt
% Mimics sntoyA.f in $SNOPT/examples
% FOR DELTA DOGS
% Example of the 'fmincon'-style call to SNOPT.
%
%     Minimize      100*x(2)
%
%     subject to             x(2)  - 1/2*(x(1)-1/2)^2 -  0.1  >= 0
%                            x(2)                             >= 0.
%
%
snscreen on;
snprint('toymin.out');  % By default, screen output is off;
sntoy.spc = which('sntoy.spc');
snspec (sntoy.spc);
snseti ('Major Iteration limit', 250);
%%
%%   control constraints
%
x0=[0.4;0];
A      = [ 0  0];
b      = [ 0 ];
Aeq    = [];
beq    = [];
%
lb     = -Inf*ones(2,1);
ub     = Inf*ones(2,1);
n= length(x0);
bnd1 = zeros(n,1);
bnd2 = ones(n,1);
inter_method=1;
% upper and lower bouds
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
%%
options.name = 'toyprob';
options.stop = @toySTOP;
[x,fval,INFO,lambda] = snsolve( @objfungrad2, x0, A, b, Aeq, beq, lb, ub, @confungrad2, options);
% [x,fval,INFO,lambda] = snsolve( @objfungrad2, x0, Ain, bin, Aeq, beq, lb, ub, @confungrad2, options);
snprint off;
snend;

function [f,gradf] = objfungrad2(x)
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
function [c,ceq,DC,DCeq] = confungrad2(x)
c(1,1) = -(x(2)-0.5*(x(1)-0.5).^2-0.1); % Inequality constraints
c(2,1) = -x(2); 
% No nonlinear equality constraints
ceq=[];
% Gradient of the constraints:
% % % if nargout > 2
    DC= [x(1)-0.5, -1;
         0, -1];
    DCeq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iAbort] = toySTOP(itn, nMajor, nMinor, condZHZ, obj, merit, step, ...
			      primalInf, dualInf, maxViol, maxViolRel, ...
			      x, xlow, xupp, xmul, xstate, ...
			      F, Flow, Fupp, Fmul, Fstate)

% Called every major iteration
% Use iAbort to stop SNOPT (if iAbort == 0, continue; else stop)

iAbort = 0


