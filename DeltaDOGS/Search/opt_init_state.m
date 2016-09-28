%initialization
function opt_init_state
% global N
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
% parameters
n=2; % design parameter dimension
ms=2; m=2*n; 
% boundary projection element
r=2;


nit = 0;
stop_opt = 0;
curr_best = [];
% for grid case
% spc = 5*ones(1,N); % initial grid points on the mesh
% delta = 1;


save opt_prev nit curr_best stop_opt 
Jall = [];
Aall = [];
save surr_pts Jall Aall

