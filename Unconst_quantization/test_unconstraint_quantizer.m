    % test the unconstraint quantizer
    % shahrouz alimo August 11 2016
    clear all;clc;
    %%
    %lattice= 'Zn ' , 'An ', 'An*?!','Dn ','Dn*', 'E8 '
    global lattice plane 
    % initialization 
    n=8;  S=rand(n,1) % data point
    %%
    lattice='E8 ';   [~,B,plane]=init_DOGS(n,lattice);
    %scale factor. Different for each lattice!
    scale=5;
    %
    %quantization 
    [Sq,Errq,Sz]=Unconstraint_quantizer(S,scale, lattice);
    
    X = [S, Sq]
    e = Errq