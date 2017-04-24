
   clear; close all;
   fprintf('implement time:%d.%d.%d.%d.%d.%d.\n',fix(clock));  % display the time
   my_import;
                                                                         
   %[syn_para] = Set_Synthetic_Data_Para;              fprintf(' finish Set_Synthetic_Data_Para \n');
   %[true_data] = Generate_Synthetic_Data(syn_para);   fprintf(' finish Generate_Synthetic_Data \n');  
   [true_data] = Import_Data;                      fprintf('\n finish Import_Data \n');
   
   [para] = Set_Algorithm_Para(20,15,1,4);             fprintf(' finish Set_Algorithm_Para \n');
   [est_data,h] = Trend_Separation(para, true_data);  fprintf('\n finish Trend_separation \n');
   save('F_data.mat','para','true_data','est_data');
   
   %Plot_Result_on_Screen(para, true_data, est_data); fprintf('\n finish Plot_Result \n');
   Output_Result_to_Figures;                          fprintf('\n finish Output_Result_to_Figures \n');
   