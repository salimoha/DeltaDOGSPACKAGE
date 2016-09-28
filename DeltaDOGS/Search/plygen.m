% before running type the following four lines to set up the environment
% setenv('DYLD_LIBRARY_PATH', '/opt/local/bin:/opt/local/lib:');                 % required for AVL to use the system Fortran library instead of the Matlab ones
% if count(py.sys.path,'') == 0                                                  % reads available python modules in the folder
% insert(py.sys.path,int32(0),'');
% end

function eff = plygen(params)

params=params.';
params=[params 0 0 0 0];
pyparams = py.list(params);

status = 1;
Nchord = 12;
Nspan = 101;
while status ~= 0 
  py.plygen.preAVL(pyparams,Nchord,Nspan);                                                   % preprocessing: writing the AVL input files
  
  cmd = './avl foil.avl < avl_script' ;
  [status,cmdout] = system(cmd);                                                 % processing: run avl and grep the output
  Nchord = Nchord + 1;
  Nspan = Nspan + 1;

  status
end
A = strsplit(cmdout,'\n');                                                     % convert the output to an array separtating the lines
cmdout = py.list(A);                                                           % convert the output to a python list

eff = py.plygen.postAVL(cmdout);
eff=-1/eff;

end
