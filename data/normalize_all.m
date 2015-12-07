clc;
clear all;
close all;

files = {'AP1','AP2','AP3','AP4','AP5','AP6','AP7','AP8','AP9'...
          'BP1','BP2','BP3','BP4','BP5','BP6','BP7','BP8','BP9','BP10'...
          'CP1','CP2','CP3','CP4','CP5','CP6','CP7','CP8','CP9','CP10','CP11'...
          'DP1','DP2','DP3','DP4','DP5','DP6','DP7','DP8','DP9','DP10','DP11'...
          'EP1','EP2','EP3','EP4','EP5','EP6','EP7','EP8','EP9','EP10','EP11'...
          'FP1','FP2','FP3','FP4','FP5','FP6','FP7','FP8'...
          'GP1','GP2','GP3','GP4','GP5','GP6','GP7','GP8','GP9','GP10','GP11'...
          'HP1','HP2','HP3','HP4','HP5','HP6','HP7','HP8','HP9','HP10','HP11'...
          'IP1','IP2','IP3','IP4','IP5','IP6','IP7','IP8','IP9','IP10','IP11'};
      
 l = length(files);
 train_data = [];
 
 for i = 1:l
     load(char(files(i)));
     
     train_data = [train_data' master_trainer']';
     
 end
 
 mn = mean(train_data);
 mx = (max(train_data));
 mx = mx + .000000001*ones(size(mx));
 
 for i = 1:l
     load(char(files(i)));
     
     s = size(master_trainer);
     master_trainer = master_trainer - ones(s(1),1)*mn;
     
     mx = abs(mx - mn.*ones(size(mx)));
     
     master_trainer = 100*master_trainer./( ones(s(1),1)*mx );
     
     save(char(files(i)),'master_trainer','responsevar','var_x','mean_x',...
                         'diff_x', 'AR_coef', ...
                         'range');
 end