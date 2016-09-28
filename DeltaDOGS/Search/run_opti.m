% function run_opti

warning off;


opt_init();
stop=0;
% while (stop == 0)
for stop =1
%     DOGS_stand_alone()
%     test_nDOGS_stand_alone()
    f=fopen('store_neval','r');
    neval=str2num(fgetl(f));
    fclose(f);
    f1=fopen('pts_to_eval','r');
    f2=fopen('surr_yi_new','w');
    f3=fopen('surr_xi_new','w');
    for i=1:neval
        rline=[];
        tmp=fgetl(f1);
        rline=[rline;str2num(tmp)]
        ym = sum(rline.^2);
        xx = str2num(tmp);
    
      fprintf(f2,'%12f\n',ym);
%       fprintf(f3,'%12f\n',rline(1));
%       fprintf(f3,'%12f %12f \n',rline(1),rline(2));
keyboard
       fprintf(f3,'%12f %12f %12f %12f\n',rline(1),rline(2),rline(3),rline(4));
%       fprintf(f3,'%12f %12f %12f %12f \n',xx(1),xx(2),xx(3),xx(4));
   end
    fclose(f1); 
    fclose(f2); 
    fclose(f3);
    f=fopen('stop_file','r');
    stop=str2num(fgetl(f));
    fclose(f);
end
