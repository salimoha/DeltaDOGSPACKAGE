function [xm cse indm ym] = tringulation_search_bound1(inter_par,xi,tri)

global n bnd1 bnd2 y0
% Update global interpolation
    %xie=xiT(:,yiT<y_cr); yie=yiT(:,yiT<y_cr);
    %[w,v] = polyharmsp_weight3(xie, yie.',1); tri=delaunayn(xiT.');
ym=inf; cse=2;
    for ind=1:size(tri,1)
      [xc,R2]=circhyp(xi(:,tri(ind,:)), n);
      if R2<5
       x=xi(:,tri(ind,:))*ones(n+1,1)/(n+1);
       y=(interpolate_val(x,inter_par)-y0)/(R2-norm(x-xc)^2);
      if (y<ym)
          ym=y; xm=x; indm=ind;
      end
      end
    end
    
    [xc,R2]=circhyp(xi(:,tri(indm,:)), n);
    x=xi(:,tri(indm,:))*ones(n+1,1)/(n+1); 
    [x y cse]=Adoptive_K_Search(x,inter_par,xc,R2);
    if cse==1
       xm=inter_min(xm,inter_par);
    end
end