function [xm ym CSm indm] = tringulation_search_constraints(inter_par_p,inter_par_g,xi)
% This function evalutes the 
global n bnd1 bnd2 Search tri
tri=delaunayn(xi.');
% keyboard
ym=inf;
CSm=1;
%            keyboard
    for ind=1:size(tri,1)
%          keyboard
      [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
      if sqrt(R2)<2e3
          %%%% ??? what does it do?
          x=xi(:,tri(ind,:))*ones(n+1,1)/(n+1);
        if Search.method ==1
%                keyboard
if length(Search.constant) ~=1

 [x y CS]=Constant_K_Search_FilterSQP_multi_K(x,inter_par_p, inter_par_g,xc,R2,Search);
else
%             [x y CS]=Constant_K_Search_FilterSQP(x,inter_par_p, inter_par_g,xc,R2,Search);
           %     keyboard
%          [x y CS]=Constant_K_Search(x,inter_par_p,inter_par_g,xc,R2,Search);
%         keyboard
%          [x y CS]=Constant_K_Search_sqp(x,inter_par_p,inter_par_g,xc,R2, Search);
%        keyboard
      [x y CS]=Constant_K_Search_snopt(x,inter_par_p,inter_par_g,xc,R2,Search);
end
         if (CS<CSm)
             xm=x; ym=y; CSm = CS; indm=ind;
         end
         if (y<ym & CS==CSm)
%          elseif (y<ym & CS==CSm)
             xm=x; ym=y; indm=ind;
         end

         %%
         
%           keyboard
          elseif Search.method ==2  
               ind
%                 keyboard
          [x y CS]=Constraint_Adaptive_K_Search(x,inter_par_p,inter_par_g,xc,R2, Search); 
%           [x y CS]=Constraint_Adaptive_K_Search_inverse(x,inter_par_p,inter_par_g,xc,R2, Search); 
%             keyboard
      if CS ==0
       % there is a point that pk < y0
           xm=x; ym=y; CSm = CS; indm=ind;
           break
      end
      if y < ym & CS ==1
       xm=x; ym=y; CSm = CS; indm=ind;      
      end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
      end   
% go to the next simplex till you find a point that satisfies the
% abovementioned criteria
      end
    end
    
