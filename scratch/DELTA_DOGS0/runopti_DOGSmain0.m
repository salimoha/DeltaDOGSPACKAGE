
fun=@(x) 5*norm((x-0.3)).^2;
yE=[];
for k=0:5
clear yE xE
    mainDOGS0(k)
    
    % load the DOGS suggestions
    xE = load('pts_to_eval')
    % function evalution 
    for jj=1:size(xE,2),  ym(jj) = fun(xE(:,jj)); end
    
     save surr_J_new ym -ASCII
    if k==0
        yE=[];
    else
    yE=load('surr_J_eval_set')
    end    
yE=[yE,ym];
    save 'surr_J_eval_set' yE -ASCII
%      end


keyboard
end
