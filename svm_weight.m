clc;
clear all;
close all;

win_time = 1;
trainer_extraction;

weight = zeros(1,9);
         
load(['train_data' '_' num2str(win_time*1000)]);
load(['test_data' '_' num2str(win_time*1000)]);

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
     
     %l = size(master_trainer);
    
     [predicted_label, accuracy, decision_values] = svmpredict(res_test, ...
                                                    test_data,model);
                      
     %result1 = length( find(predicted_label == responsevar(1)));
     %fprintf('testing data from %d\n',j);
    % fprintf('total match : %d\n',result1);
     %fprintf('total data : %d\n',l(1));
     
      result = length( find(res_test == 1));
      result = [result length( find(res_test == 2))];
      result = [result length( find(res_test == 3))];
      result = [result length( find(res_test == 4))];
      result = [result length( find(res_test == 5))];
      result = [result length( find(res_test == 6))];
      result = [result length( find(res_test == 7))];
      result = [result length( find(res_test == 8))];
      given_data = [result length( find(res_test == 9))]
     
      fprintf('**************************************\n');
           
      ln = round(sum(given_data)/9);
      
      for i =1:9
          data = predicted_label((i-1)*ln+1:i*ln);
          result = length( find(data == 1));
          result = [result length( find(data == 2))];
          result = [result length( find(data == 3))];
          result = [result length( find(data == 4))];
          result = [result length( find(data == 5))];
          result = [result length( find(data == 6))];
          result = [result length( find(data == 7))];
          result = [result length( find(data == 8))];
          result = [result length( find(data == 9))];
          
          weight(i) = result(i);
          
          fprintf('grid %d result : ',i);
          fprintf('%4d ',result);
          fprintf('\n');
      end
     
     %[~,grid] = max(result)

     fprintf('**************************************\n');
     
     weight = 100*weight./max(weight)
     
     save(['weight' '_' num2str(win_time*1000)]);
 %end
 
% fprintf('Efficiency is %f%%\n',success*100/len(2));