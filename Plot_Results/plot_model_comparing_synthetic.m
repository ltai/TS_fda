% plot model comparing, designed for synthetic data

% created on Jan.21, 2016
% continued on Feb. 16, 2016;
% continued on Feb. 21, 2016; finish all the output, save to png and txt

clear; close all;
load model_comparisom.mat;

markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};
colors = {'b','r','g','y','c','m','k'};

%%%%%%%%% plot all g %%%%%%%%
fig_all_g = figure(100); My_Figure(11,6);
for i = 1:6
    plot(true_data.t,est_g(:,i),'linewidth',1,'marker',markers{i},'color',colors{i}); hold on;
end
plot(true_data.t, true_data.g,'k','linewidth',2);
g_legend = legend('Simple Sep.','Sep. L2','Sep. Alignment',...
    'Sep. Same Warp','Simple Est.','Alignment','true g','location','eastoutside');
set(g_legend,'FontSize',18);
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');
title('\fontsize{30}comparison of g');

%%%%%%%%%%%% plot all h %%%%%%%%%%%%%%
fig_all_h = figure(101); My_Figure(7,6);
for i = 1:4
    plot(true_data.t,est_h(:,i),'linewidth',1,'marker',markers{i},'color',colors{i}); hold on;
end
plot(true_data.t, true_data.h,'k','linewidth',2);
h_legend = legend('Simple Sep.','Sep. L2','Sep. Alignment',...
    'Sep. Same Warp','true h','location','best'); 
set(h_legend,'FontSize',18);
xlabel('\fontsize{30}time');
ylabel('\fontsize{30}value');
title('\fontsize{30}comparison of h');


%%%%%%%%%%%%% plot all true data %%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot all true data %%%%%%%%%%%%%%
    
    fig_g = figure(1); My_Figure(7,6);
    plot(true_data.t,true_data.g,'linewidth',2);
    title('\fontsize{36}true g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_r = figure(2); My_Figure(5.8,6);
    plot(true_data.t,true_data.r,'linewidth',2);
    title('\fontsize{36}warping functions');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}time');
    
    fig_h = figure(3); My_Figure(7,6);
    plot(true_data.t,true_data.h,'linewidth',2)
    title('\fontsize{36}true trend h');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_f = figure(4); My_Figure(7,6);
    plot(true_data.t,true_data.f,'linewidth',2);
    title('\fontsize{36}observed signals f_{i}(t)');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_mean_f = figure(5); My_Figure(7,6);
    plot(true_data.t,mean(true_data.f,2),'linewidth',2);
    title('\fontsize{30}mean of observed signals f_{i}(t)');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');

%%%%%%%% output estimating data %%%%%%%%%%%
%%%%%%%% output estimating data %%%%%%%%%%%

% 1 g,h
    fig_1_est_g = figure(10); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,1),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_1_est_h = figure(11); My_Figure(7,6);
    plot(true_data.t,true_data.h,true_data.t,est_h(:,1),'--r','linewidth',2);
    legend('\fontsize{26}true h','\fontsize{26}est. h','location','Best');
    title('\fontsize{36}true h and est. h');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 
    
    fig_1_L2cost = figure(12); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,1)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
% 2 g,r,h
    fig_2_est_g = figure(20); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,2),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_2_est_h = figure(21); My_Figure(7,6);
    plot(true_data.t,true_data.h,true_data.t,est_h(:,2),'--r','linewidth',2);
    legend('\fontsize{26}true h','\fontsize{26}est. h','location','Best');
    title('\fontsize{36}true h and est. h');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 

    fig_2_gam = figure(22); My_Figure(5.8,6);
    plot(true_data.t,gam(:,:,2),'linewidth',2);  
    title('\fontsize{36}est. gamma');
    xlabel('\fontsize{36}time'); 
    ylabel('\fontsize{36}time'); 
    
    fig_2_L2cost = figure(23); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,2)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
% 3 g,r,h
    fig_3_est_g = figure(30); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,3),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_3_est_h = figure(31); My_Figure(7,6);
    plot(true_data.t,true_data.h,true_data.t,est_h(:,3),'--r','linewidth',2);
    legend('\fontsize{26}true h','\fontsize{26}est. h','location','Best');
    title('\fontsize{36}true h and est. h');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 

    fig_3_gam = figure(32); My_Figure(5.8,6);
    plot(true_data.t,gam(:,:,2),'linewidth',2);  
    title('\fontsize{36}est. gamma');
    xlabel('\fontsize{36}time'); 
    ylabel('\fontsize{36}time'); 

    fig_3_L2cost = figure(33); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,3)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
% 4 g,r,h
    fig_4_est_g = figure(40); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,4),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_4_est_h = figure(41); My_Figure(7,6);
    plot(true_data.t,true_data.h,true_data.t,est_h(:,4),'--r','linewidth',2);
    legend('\fontsize{26}true h','\fontsize{26}est. h','location','Best');
    title('\fontsize{36}true h and est. h');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 

    fig_4_gam = figure(42); My_Figure(5.8,6);
    plot(true_data.t,gam(:,:,4),'linewidth',2);  
    title('\fontsize{36}est. gamma');
    xlabel('\fontsize{36}time'); 
    ylabel('\fontsize{36}time'); 

    fig_4_L2cost = figure(43); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,4)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
% 5 g
    fig_5_est_g = figure(50); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,5),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
    fig_5_L2cost = figure(51); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,5)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 

% 6 g,r
    fig_6_est_g = figure(60); My_Figure(7,6);
    plot(true_data.t,true_data.g,true_data.t,est_g(:,6),'--r','linewidth',2);
    legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
    title('\fontsize{36}true g and est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');

    fig_6_gam = figure(62); My_Figure(5.8,6);
    plot(true_data.t,gam(:,:,6),'linewidth',2);  
    title('\fontsize{36}est. gamma');
    xlabel('\fontsize{36}time'); 
    ylabel('\fontsize{36}time');
    
    fig_6_L2cost = figure(63); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,6)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
root_dir = pwd;
if( isdir('html') )
    rmdir('html','s');
end
mkdir('html');
cd('html');

print(fig_all_g,'-dpng','all_g');
print(fig_all_h,'-dpng','all_h');

print(fig_f,'-dpng','f'); 
print(fig_g,'-dpng','g');
print(fig_h,'-dpng','h');
print(fig_r,'-dpng','r');
print(fig_mean_f,'-dpng','mean_f');

print(fig_1_est_g,'-dpng','1_est_g'); 
print(fig_1_est_h,'-dpng','1_est_h');       
print(fig_1_L2cost,'-dpng','1_L2cost');   

print(fig_2_est_g,'-dpng','2_est_g'); 
print(fig_2_est_h,'-dpng','2_est_h');      
print(fig_2_gam,'-dpng','2_gam'); 
print(fig_2_L2cost,'-dpng','2_L2cost');  

print(fig_3_est_g,'-dpng','3_est_g'); 
print(fig_3_est_h,'-dpng','3_est_h');      
print(fig_3_gam,'-dpng','3_gam'); 
print(fig_3_L2cost,'-dpng','3_L2cost');  

print(fig_4_est_g,'-dpng','4_est_g'); 
print(fig_4_est_h,'-dpng','4_est_h');      
print(fig_4_gam,'-dpng','4_gam'); 
print(fig_4_L2cost,'-dpng','4_L2cost');  

print(fig_5_est_g,'-dpng','5_est_g'); 
print(fig_5_L2cost,'-dpng','5_L2cost');  

print(fig_6_est_g,'-dpng','6_est_g'); 
print(fig_6_gam,'-dpng','6_gam'); 
print(fig_6_L2cost,'-dpng','6_L2cost');  

%%%%%% output data to txt files %%%%%%%
%%%%%% output data to txt files %%%%%%%
fid = fopen('Result.txt','wt');
fprintf(fid,'implement time:%d.%d.%d.%d.%d.%d\n\n',fix(clock));
fprintf(fid, 'T: %d, N: %d\n',true_data.T,true_data.N);
if para.use_data == 1
    fprintf(fid,'generate g: %d\n',para.generate_g);
    fprintf(fid,'generate h: %d\n',para.generate_h);
end
fprintf(fid,'basis: 1 is Legen; 2 is MFS; 3 is Sin; 4 is Cos\n');
fprintf(fid,'number is from 1 to max %d \n\n',max_h_index);
fprintf(fid,'true ortho: %e\n',true_ortho);

fprintf(fid, 'idx2 (SSL2) basis: %d\n',idx2_basis);
fprintf(fid, 'idx2 (SSL2) number: %d\n',idx2_number);
fprintf(fid, 'idx3 (SS-align) basis: %d\n',idx3_basis);
fprintf(fid, 'idx3 (SS-align) number: %d\n\n',idx3_number);

%%% 1: g,h
fprintf(fid,'1, simple separation\n');
fprintf(fid,'ortho: %e\n',ortho(1));
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(1));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(1));
    fprintf(fid,'h error: %e\n',L2_h_error(1));
end

%%% 2: g,r,h
fprintf(fid,'\n2, SS-L2\n');
fprintf(fid,'ortho: %e\n',ortho(2));
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(2));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(2));
    fprintf(fid,'r error: %e\n',L2_r_error(2));
    fprintf(fid,'h error: %e\n',L2_h_error(2));
end

%%% 3: g,r,h
fprintf(fid,'\n3, SS-align\n');
fprintf(fid,'ortho: %e\n',ortho(3));
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(3));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(3));
    fprintf(fid,'r error: %e\n',L2_r_error(3));
    fprintf(fid,'h error: %e\n',L2_h_error(3));
end

%%% 4: g,r,h
fprintf(fid,'\n4, SS-same warp\n');
fprintf(fid,'ortho: %e\n',ortho(4));
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(4));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(4));
    fprintf(fid,'r error: %e\n',L2_r_error(4));
    fprintf(fid,'h error: %e\n',L2_h_error(4));
end

%%% 5: g
fprintf(fid,'\n5, simple estimation\n');
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(5));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(5));
end

%%% 6: g,h
fprintf(fid,'\n6, alignment\n');
fprintf(fid,'Final L2 cost: %e\n',final_L2cost(6));
if para.use_data == 1
    fprintf(fid,'g error: %e\n',L2_g_error(6));
    fprintf(fid,'r error: %e\n',L2_r_error(6));
end
fclose(fid);

save model_comparisom.mat;

cd(root_dir);
close all;

