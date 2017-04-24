% This is a program for generating sequence of trend extraction results
%
% created at Nov. 23, 2016
clear; close all;
fprintf('implement time:%d.%d.%d.%d.%d.%d\n',fix(clock));  % display the time
my_import;

% import data and compute
max_iter = 3;

[syn_para] = Set_Synthetic_Data_Para; 
[true_data] = Generate_Synthetic_Data(syn_para);

[para1] = Set_Algorithm_Para(max_iter, 30 ,1,1); 
[est_data1] = Simple_Separation(para1, true_data); 

[para2] = Set_Algorithm_Para(max_iter, 30 ,1,2); 
[est_data2] = Simple_Separation(para2, true_data); 

[para3] = Set_Algorithm_Para(max_iter, 25 ,1,3); 
[est_data3] = Simple_Separation(para3, true_data); 

[para4] = Set_Algorithm_Para(max_iter, 20 ,1,2); 
[est_data4] = Simple_Separation(para4, true_data); 

[para5] = Set_Algorithm_Para(max_iter, 15 ,1,2); 
[est_data5] = Simple_Separation(para5, true_data); 

% plot
fig_h1 = figure(1); My_Figure(7,6);
plot(true_data.t,true_data.h,true_data.t,est_data1.est_h(:,end),'--r','linewidth',2);
legend('true h','est h','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_h2 = figure(2); My_Figure(7,6);
plot(true_data.t,true_data.h,true_data.t,est_data2.est_h(:,end),'--r','linewidth',2);
legend('true h','est h','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_h3 = figure(3); My_Figure(7,6);
plot(true_data.t,true_data.h,true_data.t,est_data3.est_h(:,end),'--r','linewidth',2);
legend('true h','est h','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');
 
fig_h4 = figure(4); My_Figure(7,6);
plot(true_data.t,true_data.h,true_data.t,est_data4.est_h(:,end),'--r','linewidth',2);
legend('true h','est h','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_h5 = figure(5); My_Figure(7,6);
plot(true_data.t,true_data.h,true_data.t,est_data5.est_h(:,end),'--r','linewidth',2);
legend('true h','est h','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');
%
fig_g1 = figure(6); My_Figure(7,6);
plot(true_data.t,true_data.g,true_data.t,est_data1.est_g(:,end),'--r','linewidth',2);
legend('true g','est g','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_g2 = figure(7); My_Figure(7,6);
plot(true_data.t,true_data.g,true_data.t,est_data2.est_g(:,end),'--r','linewidth',2);
legend('true g','est g','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_g3 = figure(8); My_Figure(7,6);
plot(true_data.t,true_data.g,true_data.t,est_data3.est_g(:,end),'--r','linewidth',2);
legend('true g','est g','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');
 
fig_g4 = figure(9); My_Figure(7,6);
plot(true_data.t,true_data.g,true_data.t,est_data4.est_g(:,end),'--r','linewidth',2);
legend('true g','est g','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

fig_g5 = figure(10); My_Figure(7,6);
plot(true_data.t,true_data.g,true_data.t,est_data5.est_g(:,end),'--r','linewidth',2);
legend('true g','est g','location','best');
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');

% output figures
root_dir = pwd;
if( isdir('html') )
    rmdir('html','s');
end
mkdir('html');
cd('html');

print(fig_h1,'-dpng','fig_h1'); 
print(fig_h2,'-dpng','fig_h2'); 
print(fig_h3,'-dpng','fig_h3'); 
print(fig_h4,'-dpng','fig_h4'); 
print(fig_h5,'-dpng','fig_h5'); 

print(fig_g1,'-dpng','fig_g1'); 
print(fig_g2,'-dpng','fig_g2'); 
print(fig_g3,'-dpng','fig_g3'); 
print(fig_g4,'-dpng','fig_g4'); 
print(fig_g5,'-dpng','fig_g5'); 

save simple_projection;
cd(root_dir);
close all;
