% this one is for artificial, model comparing, multiple iterations
% created at Jan. 20, 2016; multiple iterations on number of basis
% finished on Feb. 16, 2016; done; coed is very robust. 
%                            multiple iterations on different basis
% continue on Feb. 21, 2016; add output error of g,r,h,L2,final L2, ortho. 
% continue on July,6,2016; modified to new data structure; done;
%                          future work: can run a particular model

clear; close all;
fprintf('implement time:%d.%d.%d.%d.%d.%d\n',fix(clock));  % display the time
my_import;

%%%%%%% multiple iteration to pick up index for h   
min_h_index = 1;  % always starts from 1, if 2, go to Set_Algorithm_para
max_h_index = 6;
max_iter = 20;

for basis_type = 1:4                 % from 1 to 4 are Cos, Sin, MFS,sLegen type
    for basis_number = 1:max_h_index   % from 1,2 to 10
        
        [syn_para] = Set_Synthetic_Data_Para;             fprintf(' finish Set_Synthetic_Data_Para \n');
        [true_data] = Generate_Synthetic_Data(syn_para);  fprintf(' finish Generate_Synthetic_Data \n');  
        %[true_data] = Import_Data;                       fprintf('\n finish Import_Data \n');
        
        [para] = Set_Algorithm_Para(max_iter, 10+5*basis_type ,min_h_index,basis_number);  
        fprintf(' finish Set_Algorithm_Para \n');
        
        %[est_data1(basis_number,basis_type)] = Simple_Separation(para, true_data); 
        [est_data2(basis_number,basis_type)] = Trend_Separation(para, true_data);
        %[est_data3(basis_number,basis_type)] = Separation_Align_Warping(para, true_data);
        %[est_data4(basis_number,basis_type)] = Separation_Same_Warping(para, true_data);
        
        %[est_data3(k2,k1)] = Separation_Area_Invariance(para, true_data);  no need? 
    end %%%%%%%%%% end of multiple iteration of h index
end %%%%%%%%%% end of multiple h basis
    %[est_data5] = Simple_Estimation(para, true_data);
    %[est_data6] = Alignment_Model(para, true_data);

%%% now is pick up the index for h %%%    
    est_g = zeros(true_data.T,6);
    est_h = zeros(true_data.T,4);
    gam = zeros(true_data.T,true_data.N,6);
    L2cost = zeros(para.max_iter,6);
    
    final_L2cost = zeros(1,6);
    L2_g_error = zeros(1,6);
    L2_h_error = zeros(1,6);
    L2_r_error = zeros(1,6);
    ortho = zeros(1,6);
    
    %temp1 = cell2mat({est_data1.final_L2cost}); 
    temp2 = cell2mat({est_data2.final_L2cost}); 
    [~,idx2] = min(temp2);
    idx2_basis = ceil(idx2/max_h_index);
    idx2_number = idx2 - max_h_index*(idx2_basis - 1);
    
   % temp3 = cell2mat({est_data3.final_SRVF_L2cost}); 
   % [~,idx3] = min(temp3);
   % idx3_basis = ceil(idx3/max_h_index);
   % idx3_number = idx3 - max_h_index*(idx3_basis - 1);
   % temp4 = cell2mat({est_data4.final_L2cost}); 
      
   % est_g(:,1) = est_data1(idx2_number,idx2_basis).est_g(:,end); % borrow L2-model
    est_g(:,2) = est_data2(idx2_number,idx2_basis).est_g(:,end);
   % est_g(:,3) = est_data3(idx3_number,idx3_basis).est_g(:,end);
   % est_g(:,4) = est_data4(idx2_number,idx2_basis).est_g(:,end); % borrow L2-model
   % est_g(:,5) = est_data5.est_g(:,end);
   % est_g(:,6) = est_data6.est_g(:,end);
    
   % est_h(:,1) = est_data1(idx2_number,idx2_basis).est_h(:,end);
    est_h(:,2) = est_data2(idx2_number,idx2_basis).est_h(:,end);
   % est_h(:,3) = est_data3(idx3_number,idx3_basis).est_h(:,end);
   % est_h(:,4) = est_data4(idx2_number,idx2_basis).est_h(:,end);
    
    gam(:,:,2) = est_data2(idx2_number,idx2_basis).gam(:,:,end);
   % gam(:,:,3) = est_data3(idx3_number,idx3_basis).gam(:,:,end);
   % gam(:,:,4) = est_data4(idx2_number,idx2_basis).gam(:,:,end);
   % gam(:,:,6) = est_data6.gam(:,:,end);
    
   % L2cost(:,1) = est_data1(idx2_number,idx2_basis).L2cost;
    L2cost(:,2) = est_data2(idx2_number,idx2_basis).L2cost(3,:)';
   % L2cost(:,3) = est_data3(idx3_number,idx3_basis).SRVF_L2cost;
   % L2cost(:,4) = est_data4(idx2_number,idx2_basis).SRVF_L2cost;
   % L2cost(:,5) = est_data5.L2cost;
   % L2cost(:,6) = est_data6.SRVF_L2cost;
    
   % final_L2cost(1) = L2cost(end,1);
    final_L2cost(2) = L2cost(end,2);
   % final_L2cost(3) = est_data3(idx3_number,idx3_basis).L2cost(end);
   % final_L2cost(4) = est_data4(idx2_number,idx2_basis).L2cost(end);
   % final_L2cost(5) = L2cost(end,5);
   % final_L2cost(6) = est_data6.L2cost(end);
    
   % ortho(1) = est_data1(idx2_number,idx2_basis).ortho;
    ortho(2) = est_data2(idx2_number,idx2_basis).ortho(end);
   % ortho(3) = est_data3(idx3_number,idx3_basis).ortho;
   % ortho(4) = est_data4(idx2_number,idx2_basis).ortho;
    
    if true_data.use_data == 1
        true_ortho = Check_Orthogonality(true_data.g,true_data.h,true_data.t);
        
       % L2_g_error(1) = est_data1(idx2_number,idx2_basis).g_error;
        L2_g_error(2) = est_data2(idx2_number,idx2_basis).g_error(end);
       % L2_g_error(3) = est_data3(idx3_number,idx3_basis).g_error;
       % L2_g_error(4) = est_data4(idx2_number,idx2_basis).g_error;
       % L2_g_error(5) = est_data5.g_error;
       % L2_g_error(6) = est_data6.g_error;
        
       % L2_h_error(1) = est_data1(idx2_number,idx2_basis).h_error;
        L2_h_error(2) = est_data2(idx2_number,idx2_basis).h_error(end);
       % L2_h_error(3) = est_data3(idx3_number,idx3_basis).h_error;
       % L2_h_error(4) = est_data4(idx2_number,idx2_basis).h_error;
        
        L2_r_error(2) = est_data2(idx2_number,idx2_basis).r_error(end);
       % L2_r_error(3) = est_data3(idx3_number,idx3_basis).r_error;
       % L2_r_error(4) = est_data4(idx2_number,idx2_basis).r_error;
       % L2_r_error(6) = est_data6.r_error;
    end
        
    save model_comparisom.mat;
    