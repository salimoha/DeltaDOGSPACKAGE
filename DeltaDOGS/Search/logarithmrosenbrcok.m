% check logarithm Idea for the Rosenbrockn
clear all
close all
fun=@rosenbrockwithgrad;
xv=-2:0.1:2;

for ii=1:length(xv)
    for jj=1:length(xv)
    %Ul(ii,jj)=log(1+fun([xv(ii);xv(jj)]);
     U(ii,jj)=fun([xv(ii);xv(jj)]);
     if U(ii,jj)<1
         Ul(ii,jj)=U(ii,jj);
     else
        Ul(ii,jj)=2- 1/U(ii,jj);
     end
    end
end
figure(1)
subplot(1,2,1)
surf(U/max(max(U)))
view(0 ,90)
colorbar
subplot(1,2,2)

surf(Ul/max(max(Ul)))
colorbar
view(0 ,90)