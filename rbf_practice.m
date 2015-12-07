clc;
clear all;
close all;

win_time = 2;
trainer_extraction;

%trainer_extraction;

files = {'PP1','PP2','PP3','PP4','PP6','PP7','PP8','PP9',...
          'PP10','PP11','PP12','PP13','PP14','PP15','PP16','PP17','PP18','PP19',...
          'PP20','PP21','PP22','PP23','PP24','PP25','PP26','PP27','PP28','PP29'...
          'PP30','PP31','PP32','PP33','PP34','PP35','PP36','PP37','PP38','PP39'...
          'PP40','PP41','PP42','PP43','PP44','PP45','PP46','PP47','PP48','PP49','PP50'...
          };
         

 bestcv = 0;
 
 for log2c = -1.1:3.1
     for log2g = -4.1:1.1
         cmd = ['-t 2 -v 2 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g) ];
         cv = svmtrain(res_train, train_data, cmd);
         
         if cv >= bestcv
             bestcv = cv;
             bestc = 2^log2c;
             bestg = 2^log2g;
         end
     end
 end

cmd = ['-t 2 -c ',num2str(bestc), ' -g ',num2str(bestg)];
%model = svmtrain(res_train, train_data, '-t 2 -c 7.4643 -g 1.8661');
model = svmtrain(res_train, train_data, cmd);

 %len = size(files);
 %success = 0;
 

% bestc
% bestg
 fprintf('\n\n#############\n');
 %for j = 1:len(2)
     %load(char(files(j)));
     
  fprintf('**************************************\n');
    
  len = length(files);
  GRD = [];
  for j =1:len   
     load([char(files(j)) '_' '2000']);
     l = size(master_trainer);
     
     [predicted_label, accuracy, decision_values] = svmpredict(zeros(l(1),1), ...
                                                    master_trainer,model);
     
           
      ln = size(master_trainer);
      ln = ln(1);
      
          result = length( find(predicted_label == 1));
          result = [result length( find(predicted_label == 2))];
          result = [result length( find(predicted_label == 3))];
          result = [result length( find(predicted_label == 4))];
          result = [result length( find(predicted_label == 5))];
          result = [result length( find(predicted_label == 6))];
          result = [result length( find(predicted_label == 7))];
          result = [result length( find(predicted_label == 8))];
          result = [result length( find(predicted_label == 9))];
          fprintf('grid %d result : ',j);
          fprintf('%4d ',result);
          fprintf('\n');
          [~,g] = max(result);
          GRD = [GRD g+'A'-1];
          
     
          [~,grid] = max(result);
          fprintf('predicted grid: %d\n',g);

  end
 fprintf('predicted grid: %s\n',GRD);
 fprintf('**************************************\n');
% fprintf('Efficiency is %f%%\n',success*100/len(2));
%AHCFBBIBDADBDIIIAAEHBBADCEEGBDDCGGEAIHIFGEIEFAFGC