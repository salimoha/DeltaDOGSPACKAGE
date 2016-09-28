function opt_init()
global n m ms bnd1 bnd2  Search r
% read the configuration file
conf = readConf('./DOGS.conf',0);
n = str2num(conf.NumParam);
ms = str2num(conf.NumConstraints);
m=2*n;
bnd1 = str2num(conf.amin)';
bnd2 = str2num(conf.amax)';
r = str2num(conf.ConstProjc);
Search.method = str2num(conf.method);
Search.constant = str2num(conf.optVar);

    
end

