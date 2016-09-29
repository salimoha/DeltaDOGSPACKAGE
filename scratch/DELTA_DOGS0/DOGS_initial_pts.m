function [xU,xE,Ain, bin,length_scale, deltaU] = DOGS_initial_pts(bnd1,bnd2,n, delta0)
% DOGS LAMBDA PACKAGE 
% author: Shahrouz Alimohammadi
% 08/29/2016
% This function generates the initial points to start the Delta DOGS algorithm 

% related to the xE set
if nargin < 4,   delta0=0.15; end
% range of the space
length_scale=max(bnd2-bnd1); delta0=delta0*length_scale;

% Input the equality constraints
Ain=[eye(n);-eye(n)]; bin=[bnd2 ;-bnd1];

%% Calculate the Initial support trinagulation points (body center as the
% evaluaition initial point )
xU=bounds(bnd1,bnd2, n); % xU=vertex_find(Ain,bin,[],[]);
% xU=vertex_find(Ain,bin,bnd1,bnd2);
%% Calculate initial  evaluation points with delta0 distance
xE=(bnd1+0.5*(bnd2-bnd1)); % xE=(sum(xU')/(n+1))';
% initial points: the midpoint point and its neighber
for ii=1:n; e=zeros(n,1); e(ii)=1;xE(:,ii+1)=xE(:,1)+delta0*e; end

% Calculate initial unevaluated delta (distance from xE)
for ii=1:size(xU,2),  deltaU(ii)=mindis(xU(:,ii),xE); end


end