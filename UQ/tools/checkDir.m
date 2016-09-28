function checkDir(inputDir , cleanDir )
% checks the directory 
%  if it exist : do clean or not
%  if it does not exist makdir
% Shahrouz Alima oct 14 2014


if(nargin<2), cleanDir = false; end

if cleanDir == true
    
    
    doCln = input('do you want to clean the exisiting directory ? \n','s')
    
    switch doCln
        case 'yes'
        
    rmdir(inputDir,'s')
    shwmsg = [inputDir ' has been cleaned'];
    disp(shwmsg)
    mkdir(inputDir)
    
        case 'no'            
       cleanDir = false;  
    end
end


if cleanDir == false
    if(exist(inputDir,'dir'))
    error(' cannot over write a file .... ')
    end
mkdir(inputDir)
end

if cleanDir == 3
   disp('do the simulation \n')
    
end

end

