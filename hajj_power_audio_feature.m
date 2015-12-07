%% power extraction

%% merge files
clc;
clear all;
close all;
A=9;
B=10;
C=11;
D=11;
E=11;
D=11;
E=11;
F=8;
G=11;
H=11;
I=11;
%A=2;B=2;C=2;D=2;E=2;F=2;F=2;G=2;H=2;I=2;
grids=[A B C D E F G H I];
grid_name = ['A','B','C','D','E','F','G','H','I'];

init=0;

for g=1:9
    sb=grids(g);
    sp=1;
    y=0;
    for  i= 1:sb
        init=init+1;
        filename = sprintf('ENF new/P%d.mat',init)
        load(filename);
       
    fn = sprintf('hajj_train/%sP%d',grid_name(g),sp);
    %F = ENF;
    save(fn,'ENF');
    sp = sp + 1;
    end
    
end


%% extract power enf

grids = ['A','B','C','D','E','F','G','H','I'];
count  = 12;

win_time = 2; % 2 sec window time

for Grid = grids   
    for i = 1:count
       
        filename = sprintf('hajj_train/%sP%d.mat',Grid,i);
   
        if exist(filename,'file')==2
            fprintf('loading from %s\n',filename);
            file_to_save = sprintf('features_hajj/%sP%d.mat',Grid,i);
            load(filename);
            %enf = power_enf(filename,win_time);
            feature_extract(ENF,char(Grid),file_to_save);    
        else
           fprintf('%s does not exist\n',filename); 
        end 
    end  
end


%% audio extraction
grids = ['A','B','C','D','E','F','G','H','I'];

counter = 1;
for Grid = grids
   
    load('ENF new/ENF_audio.mat');
    enf = [ ENF(counter,:) ENF(counter+1,:) ];
    
    file_to_save = sprintf('features_hajj/%sA.mat',Grid);
    
     fprintf('%d ) loading audio from Grid %s\n',counter,Grid);

    feature_extract(enf,char(Grid),file_to_save);
    counter = counter + 2;
    
end




