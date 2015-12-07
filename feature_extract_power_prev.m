grids = ['A','B','C','D','E','F','G','H','I'];
count  = 12;

win_time = 2; % 2 sec window time

for Grid = grids   
    for i = 1:count
       
        filename = sprintf('/Users/mac/Desktop/sp_cup/Grid_%s/Power_recordings/Train_Grid_%s_P%d.wav',...
            Grid,Grid,i);
        
        if exist(filename,'file')==2
            fprintf('loading from %s\n',filename);
            file_to_save = sprintf('power_prev/%sP%d',Grid,i);
            enf = power_enf(filename,win_time);
            feature_extract(enf,char(Grid),file_to_save);    
        else
           fprintf('%s does not exist\n',filename); 
        end 
    end  
end