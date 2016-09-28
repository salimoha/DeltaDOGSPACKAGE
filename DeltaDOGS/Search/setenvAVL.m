% before running type the following four lines to set up the environment
%pyversion [getenv('HOME') '/anaconda/bin/python']
setenv('DYLD_LIBRARY_PATH', '/opt/local/bin:/opt/local/lib:');                 % required for AVL to use the system Fortran library instead of the Matlab ones
if count(py.sys.path,'') == 0                                                  % reads available python modules in the folder
insert(py.sys.path,int32(0),'');
end


