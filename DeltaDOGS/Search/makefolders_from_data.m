clear all
close all
clc

load 3st_new 

mainfolder=pwd;
for index=n+2:length(yiT)
a=num2str(1e8+index); mkdir(a);
    y1(index)=plygen_n(ss0+rk.*xiT(:,index));
    copyfile('foil.avl',[mainfolder '/' a]);
end