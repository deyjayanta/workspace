%load('train_data');
% clc;
% clear all;
% close all;
% 
% win_time = 2;
extract_train_data;

s = size(res_train);
n = 1 : s(1);

grid1 = train_data(res_train(n) == 1,:);
grd1_res = res_train(res_train(n) == 1,:);

[grid1, grd1_res] = random_data(grid1,grd1_res);

grid2 = train_data(res_train(n) == 2,:);
grd2_res = res_train(res_train(n) == 2,:);

[grid2, grd2_res] = random_data(grid2,grd2_res);


grid3 = train_data(res_train(n) == 3,:);
grd3_res = res_train(res_train(n) == 3,:);

[grid3, grd3_res] = random_data(grid3,grd3_res);

grid4 = train_data(res_train(n) == 4,:);
grd4_res = res_train(res_train(n) == 4,:);

[grid4, grd4_res] = random_data(grid4,grd4_res);

grid5 = train_data(res_train(n) == 5,:);
grd5_res = res_train(res_train(n) == 5,:);

[grid5, grd5_res] = random_data(grid5,grd5_res);

grid6 = train_data(res_train(n) == 6,:);
grd6_res = res_train(res_train(n) == 6,:);

[grid6, grd6_res] = random_data(grid6,grd6_res);

grid7 = train_data(res_train(n) == 7,:);
grd7_res = res_train(res_train(n) == 7,:);

[grid7, grd7_res] = random_data(grid7,grd7_res);

grid8 = train_data(res_train(n) == 8,:);
grd8_res = res_train(res_train(n) == 8,:);

[grid8, grd8_res] = random_data(grid8,grd8_res);

grid9 = train_data(res_train(n) == 9,:);
grd9_res = res_train(res_train(n) == 9,:);

[grid9, grd9_res] = random_data(grid9,grd9_res);

%len = size(grid1);
%l = floor(mn*.65);
grd1_train = grid1(1:mn,:);
%grd1_test = grid1(l:len,:);
grd1_res_tr = grd1_res(1:mn,:);
%grd1_res_tst = grd1_res(l:len,:);

%len = size(grid2);
%l = floor(mn*.65);
grd2_train = grid2(1:mn,:);
%grd2_test = grid2(l:len,:);
grd2_res_tr = grd2_res(1:mn,:);
%grd2_res_tst = grd2_res(l:len,:);

%len = size(grid3);
%l = floor(mn*.65);
grd3_train = grid3(1:mn,:);
%grd3_test = grid3(l:len,:);
grd3_res_tr = grd3_res(1:mn,:);
%grd3_res_tst = grd3_res(l:len,:);

%len = size(grid4);
%l = floor(mn*.65);
grd4_train = grid4(1:mn,:);
%grd4_test = grid4(l:len,:);
grd4_res_tr = grd4_res(1:mn,:);
%grd4_res_tst = grd4_res(l:len,:);

%len = size(grid5);
%l = floor(mn*.65);
grd5_train = grid5(1:mn,:);
%grd5_test = grid5(l:len,:);
grd5_res_tr = grd5_res(1:mn,:);
%grd5_res_tst = grd5_res(l:len,:);

%len = size(grid6);
%l = floor(mn*.65);
grd6_train = grid6(1:mn,:);
%grd6_test = grid6(l:len,:);
grd6_res_tr = grd6_res(1:mn,:);
%grd6_res_tst = grd6_res(l:len,:);

%len = size(grid7);
%l = floor(mn*.65);
grd7_train = grid7(1:mn,:);
%grd7_test = grid7(l:len,:);
grd7_res_tr = grd7_res(1:mn,:);
%grd7_res_tst = grd7_res(l:len,:);

%len = size(grid8);
%l = floor(mn*.65);
grd8_train = grid8(1:mn,:);
%grd8_test = grid8(l:len,:);
grd8_res_tr = grd8_res(1:mn,:);
%grd8_res_tst = grd8_res(l:len,:);

%len = size(grid9);
%l = floor(mn*.65);
grd9_train = grid9(1:mn,:);
%grd9_test = grid9(l:len,:);
grd9_res_tr = grd9_res(1:mn,:);
%grd9_res_tst = grd9_res(l:len,:);

train_data = [grd1_train' grd2_train' grd3_train'...
              grd4_train' grd5_train' grd6_train'...
              grd7_train' grd8_train' grd9_train']';
          
res_train = [grd1_res_tr' grd2_res_tr' grd3_res_tr'...
             grd4_res_tr' grd5_res_tr' grd6_res_tr'...
             grd7_res_tr' grd8_res_tr' grd9_res_tr']';
         
% test_data = [grd1_test' grd2_test' grd3_test'...
%               grd4_test' grd5_test' grd6_test'...
%               grd7_test' grd8_test' grd9_test']';
%           
% res_test = [grd1_res_tst' grd2_res_tst' grd3_res_tst'...
%              grd4_res_tst' grd5_res_tst' grd6_res_tst'...
%              grd7_res_tst' grd8_res_tst' grd9_res_tst']';

%  weight_train = [grd1_test' grd2_test' grd3_test'...
%               grd4_test' grd5_test' grd6_test'...
%               grd7_test' grd8_test' grd9_test']';
          
%  weight_res = [grd1_res_tst' grd2_res_tst' grd3_res_tst'...
%              grd4_res_tst' grd5_res_tst' grd6_res_tst'...
%              grd7_res_tst' grd8_res_tst' grd9_res_tst']';
         

save(['train_data' '_' num2str(win_time*1000)],'train_data','res_train');
%save(['test_data' '_' num2str(win_time*1000)],'test_data','res_test');
% save('weight_data_400','weight_train','weight_res');