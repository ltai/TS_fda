% plot model comparing, designed for real data, 
% labels will be changed from data to data

% Copy from plot_model_comparing_synthetic on Feb. 22, 2016;
% for growth velocity data! already change case by case.

clear; close all;
load model_comparisom.mat;

markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};
colors = {'b','r','g','y','c','m','k'};

true_data.t = linspace(0,18,true_data.T);

%%%%%%%%% plot all g %%%%%%%%
fig_all_g = figure(100); My_Figure(11,6);
for i = 1:6
    plot(true_data.t,est_g(:,i),'linewidth',1,'marker',markers{i},'color',colors{i}); hold on;
end
g_legend = legend('Simple Sep.','Sep. L2','Sep. Alignment',...
    'Sep. Same Warp','Simple Est.','Alignment','location','eastoutside');
set(g_legend,'FontSize',18);
xlabel('\fontsize{30}age(years)'); xlim([0 18]);
ylabel('\fontsize{30}growth velocity (cm/year)');
title('\fontsize{30}comparison of g');

%%%%%%%%%%%% plot all h %%%%%%%%%%%%%%
fig_all_h = figure(101); My_Figure(7,6);
for i = 1:4
    plot(true_data.t,est_h(:,i),'linewidth',1,'marker',markers{i},'color',colors{i}); hold on;
end
h_legend = legend('Simple Sep.','Sep. L2','Sep. Alignment',...
    'Sep. Same Warp','location','best'); 
set(h_legend,'FontSize',18);
xlabel('\fontsize{30}age(years)'); xlim([0 18]);
ylabel('\fontsize{30}growth velocity (cm/year)');
title('\fontsize{30}comparison of h');


%%%%%%%%%%%%% plot all true data %%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot all true data %%%%%%%%%%%%%%
    
    fig_f = figure(4); My_Figure(7,6);
    plot(true_data.t,true_data.f,'linewidth',2);
    title('\fontsize{30}Berkeley Male Growth Velocity');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_mean_f = figure(5); My_Figure(7,6);
    plot(true_data.t,mean(true_data.f,2),'linewidth',2);
    title('\fontsize{24}mean of Berkeley Male Growth Velocity');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');

%%%%%%%% output estimating data %%%%%%%%%%%
%%%%%%%% output estimating data %%%%%%%%%%%

% 1 g,h
    fig_1_est_g = figure(10); My_Figure(7,6);
    plot(true_data.t,est_g(:,1),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_1_est_h = figure(11); My_Figure(7,6);
    plot(true_data.t,est_h(:,1),'--r','linewidth',2);
    title('\fontsize{36}est. h');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)'); 
    
    fig_1_L2cost = figure(12); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,1)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 
    
% 2 g,r,h
    fig_2_est_g = figure(20); My_Figure(7,6);
    plot(true_data.t,est_g(:,2),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_2_est_h = figure(21); My_Figure(7,6);
    plot(true_data.t,est_h(:,2),'--r','linewidth',2);
    title('\fontsize{36}est. h');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)'); 

    fig_2_gam = figure(22); My_Figure(5.8,6);
    plot(linspace(0,1,true_data.T),gam(:,:,2),'linewidth',2);  
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
    plot(true_data.t,est_g(:,3),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_3_est_h = figure(31); My_Figure(7,6);
    plot(true_data.t,est_h(:,3),'--r','linewidth',2);
    title('\fontsize{36}est. h');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)'); 

    fig_3_gam = figure(32); My_Figure(5.8,6);
    plot(linspace(0,1,true_data.T),gam(:,:,2),'linewidth',2);  
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
    plot(true_data.t,est_g(:,4),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_4_est_h = figure(41); My_Figure(7,6);
    plot(true_data.t,est_h(:,4),'--r','linewidth',2);
    title('\fontsize{36}est. h');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)'); 

    fig_4_gam = figure(42); My_Figure(5.8,6);
    plot(linspace(0,1,true_data.T),gam(:,:,4),'linewidth',2);  
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
    plot(true_data.t,est_g(:,5),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');
    
    fig_5_L2cost = figure(51); My_Figure(7,6);
    plot((1:para.max_iter),log10(L2cost(:,5)),'-r*','linewidth',2);
    title('\fontsize{36}cost function');
    xlabel('\fontsize{36}number of iterations');
    ylabel('\fontsize{36}log square of L2 norm'); 

% 6 g,r
    fig_6_est_g = figure(60); My_Figure(7,6);
    plot(true_data.t,est_g(:,6),'--r','linewidth',2);
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)'); xlim([0 18]);
    ylabel('\fontsize{36}growth velocity (cm/year)');

    fig_6_gam = figure(62); My_Figure(5.8,6);
    plot(linspace(0,1,true_data.T),gam(:,:,6),'linewidth',2);  
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

if para.use_data == 0
    fprintf(fid,'real data selection: %d\n', para.real_data);    
end
fprintf(fid,'basis: 1 is Legen; 2 is MFS; 3 is Sin; 4 is Cos\n');
fprintf(fid,'min h index is %d \n',min_h_index);
fprintf(fid,'max h index is %d \n\n',max_h_index);

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

