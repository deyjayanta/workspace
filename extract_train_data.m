%clc;
%clear all;
%close all;

files = {'AP1','AP2','AP3','AP4','AP5','AP6','AP7','AP8','AP9',...
           'BP1','BP2','BP3','BP4','BP5','BP6','BP7','BP8','BP9','BP10',...
           'CP1','CP2','CP3','CP4','CP5','CP6','CP7','CP8','CP9','CP10','CP11'...
           'DP1','DP2','DP3','DP4','DP5','DP6','DP7','DP8','DP9','DP10','DP11'...
           'EP1','EP2','EP3','EP4','EP5','EP6','EP7','EP8','EP9','EP10','EP11'...
           'FP1','FP2','FP3','FP4','FP5','FP6','FP7','FP8'...
           'GP1','GP2','GP3','GP4','GP5','GP6','GP7','GP8','GP9','GP10','GP11'...
           'HP1','HP2','HP3','HP4','HP5','HP6','HP7','HP8','HP9','HP10','HP11'...
           'IP1','IP2','IP3','IP4','IP5','IP6','IP7','IP8','IP9','IP10','IP11'...
           'AA1','AA2','BA1','BA2','CA1','CA2','DA1','DA2','EA1','EA2',...
           'FA1','FA2','GA1','GA2','HA1','HA2','IA1','IA2'...
        };
      
 datasample(files,length(files),'Replace',false);
      
      
 l = length(files);
 train_data = [];
 res_train = [];
 ln = 0;
 vr = 1;
 
%win_time = 1;
 
 for i = 1:l
     load([char(files(i)) '_' num2str(win_time*1000)]);
     
%      if responsevar(1,1) ~= vr
%          if ln < mn
%              mn = ln;
%          end
%          
%          k = size(responsevar);
%          ln = k(1);
%      else
%          k = size(responsevar);
%          ln = ln + k(1);
%      end
%      
%      vr = responsevar(1,1);
 
 
     train_data = [train_data master_trainer'];
     res_train = [res_train responsevar'];
     
 end
 
 mn = sum(res_train == 1);
 
 for i = 2:9
     mn = min(mn,sum(res_train == i));
 end
 
 train_data = train_data';
 res_train = res_train';
 
 %save('train_data','train_data','res_train','mn');