function fun_eval()

xE =  load('pts_to_eval');
for ii=1:size(xE,2)
yE(ii)=norm(xE(:,ii));
end
save surr_J_new yE -ASCII
   
 
   %% constraint violation
   ms = 2;
con{1}=@(x) sum(3*x)+sum(1./(3*x)); 

con{2}=@(x)-( sum(3*x)+sum(1./(3*x)) ); 
for jj=1:size(xE,2)
  for ii=1:ms, 
% C(ii,jj) = con{ii}(xE(:,jj));
  
C{ii}(jj) = con{ii}(xE(:,jj));
  end
end
% save surr_C_new C -ASCII
save surr_C_new.mat C 
end