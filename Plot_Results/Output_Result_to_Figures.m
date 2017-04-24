% this is the main function for trend separation, 
% continued on old version which is created on July, 27, 2014
% created on March, 21, 2015
% est_energy for orthogonal penalty to be ajusted, July, 23, 2015
% may use Simpson's rule to replace evaluating a quadrature!!
%            current version is using spline, may to too expansive
%                                            Aug. 12, 2015
% continued on Sep. 7, 2015; add change directory and save data, pictures
% continued on Sep. 28, 2015; revise all the figure output
% ccontinued on July, 2, 2016; adapt to new data structure; successful                          

function Output_Result_to_Figures
    load F_data; 
    %my_import;
    
    [T,~] = size(true_data.f);
    max_iter = para.max_iter;
    f = true_data.f;
    set(0,'DefaultFigureVisible','off');
    warning('off', 'all');   %%% surpress warning    

%%%%%%%%%%% data generation or observation %%%%%%%
%%%%%%%%%%% data generation or observation %%%%%%%
    if true_data.use_data == 1  %%% artificial        
        t = true_data.t;
        fig_g = figure(1); My_Figure(7,6);
        plot(t,true_data.g,'linewidth',2);
        %title('\fontsize{30}true g');
        xlabel('\fontsize{30}time'); 
        ylabel('\fontsize{30}value');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        
        fig_r = figure(2); My_Figure(5.8,6);
        plot(t,true_data.r,'linewidth',2);
        %title('\fontsize{30}warping functions');
        xlabel('\fontsize{30}time'); 
        ylabel('\fontsize{30}time');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        
        fig_h = figure(3); My_Figure(7,6);
        plot(t,true_data.h,'linewidth',2);
        %title('\fontsize{30}true trend h');
        xlabel('\fontsize{30}time');
        ylabel('\fontsize{30}value');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        
        fig_f = figure(4); My_Figure(7,6);
        plot(t,true_data.f,'linewidth',2);
        %title('\fontsize{30}observed signals f_{i}(t)');
        xlabel('\fontsize{30}time');
        ylabel('\fontsize{30}value');  
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
   
    elseif true_data.use_data == 0  %%% real data
        t = linspace(0,11,T)';  %%%% for Berkely growth velocity
        
        fig_f = figure(4); My_Figure(7,6); 
        plot(t,true_data.f,'linewidth',2); 
        xlim([0 11]); ylim([min(min(f)) max(max(f))]);
        %title('\fontsize{30}male growth velocity f_{i}(t)'); 
        xlabel('\fontsize{30}hours');
        ylabel('\fontsize{30}gene expression level'); 
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    end
    
%%%%%%%%%%%   plot estimated solutions   %%%%%%%%%%%%%%%  
%%%%%%%%%%%   plot estimated solutions   %%%%%%%%%%%%%%%   
    est_g = est_data.est_g;
    est_h = est_data.est_h;
    gam = est_data.gam;   

    fig_est_g = figure(5); My_Figure(7,6);    
    if true_data.use_data == 0
        plot(t,est_g(:,max_iter),'--r','linewidth',2); 
        xlim([0 11]); ylim([min(est_g(:,end)) max(est_g(:,end))]);
        %title('\fontsize{30}est. g');
        xlabel('\fontsize{30}hours');
        ylabel('\fontsize{30}gene expression level'); 
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    elseif true_data.use_data == 1
        plot(t,true_data.g,t,est_g(:,max_iter),'--r','linewidth',2);
        legend('\fontsize{26}true g','\fontsize{26}est. g','location','Best');
        %title('\fontsize{30}true g and est. g');
        xlabel('\fontsize{30}time'); 
        ylabel('\fontsize{30}value');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    end

    fig_gam = figure(6); My_Figure(5.8,6);
    plot(linspace(0,1,T)',gam(:,:,max_iter),'linewidth',2);  
    %title('\fontsize{30}est. gamma');
    xlabel('\fontsize{30}time'); 
    ylabel('\fontsize{30}time'); 
    ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    
    fig_est_h = figure(7); My_Figure(7,6);
    if true_data.use_data == 0
        plot(t,est_h(:,max_iter),'--r','linewidth',2);  
        xlim([0 11]); ylim([min(min(f)) max(max(f))]);
        %title('\fontsize{30}est. h');
        xlabel('\fontsize{30}hours');
        ylabel('\fontsize{30}gene expression level');  
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    elseif true_data.use_data == 1
        plot(t,true_data.h,t,est_h(:,max_iter),'--r','linewidth',2);  
        legend('\fontsize{26}true h','\fontsize{26}est. h','location','Best');
        %title('\fontsize{30}true h and est. h');
        xlabel('\fontsize{30}time'); 
        ylabel('\fontsize{30}value');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    end
   
%%%%%%%%%   plot error   %%%%%%%%%%    
%%%%%%%%%   plot error   %%%%%%%%%%
    fig_L2cost = figure(8); My_Figure(7,6);
    plot((1:max_iter), (est_data.L2cost(3,:)),'r','linewidth',2); hold on;
    plot((1:max_iter), (est_data.L2cost(3,:)),'-r*','linewidth',2);
    xlabel('\fontsize{30}number of iterations');
    ylabel('\fontsize{30}L^2 cost'); 
    %title('\fontsize{30}minimized L^2 cost');
    ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    
    %str = sprintf('L2 cost: %0.2e',est_data.L2cost(3,end));
    %temp_text = text(max_iter/3, 0.5*min(log10(est_data.L2cost(3,1:end))),str);
    %set(temp_text,'fontsize',32);

    if true_data.use_data == 1 % synthetic   
        fprintf('\ntrue sum square L2norm %e\n',est_data.syn_L2cost);
        fprintf('est. sum square L2norm %e\n',est_data.final_L2cost);
        fprintf('abslote L2_error of r %e\n', est_data.r_error(max_iter));
        fprintf('relative L2_error of g %e\n', est_data.g_error(max_iter) );
        fprintf('relative L2_error of h %e\n\n', est_data.h_error(max_iter) );
        fprintf('true ortho. %e\n', est_data.true_ortho );
        fprintf('est ortho. %e\n', est_data.ortho(max_iter));

        fig_L2cost_gam = figure(9); My_Figure(7,6);
        plot((1:max_iter), (est_data.r_error),'r','linewidth',2);
        plot((1:max_iter), (est_data.r_error),'-r*','linewidth',2);
        %title('\fontsize{30}sum of L^2 error in r_i');
        xlabel('\fontsize{30}number of iterations');
        ylabel('\fontsize{30}L2 error');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        %str = sprintf('L^2 error: %0.2e',est_data.r_error(max_iter));
        %temp_text = text(max_iter/3, 0.5*min(log10(est_data.r_error(2:end))),str);
        %set(temp_text,'fontsize',32);
        
        fig_L2cost_est_g = figure(10); My_Figure(7,6);
        plot((1:max_iter), (est_data.g_error),'r','linewidth',2 );
        plot((1:max_iter), (est_data.g_error),'-r*','linewidth',2 );
        %title('\fontsize{30}L^2 error in g');
        xlabel('\fontsize{30}number of iterations');
        ylabel('\fontsize{30}log of L2 error');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        %str = sprintf('L^2 error: %0.2e',est_data.g_error(max_iter));
        %temp_text = text(max_iter/3, 0.5*min(log10(est_data.g_error(2:end))),str);
        %set(temp_text,'fontsize',32);

        fig_L2cost_est_h = figure(11); My_Figure(7,6);
        plot((1:max_iter), (est_data.h_error),'r','linewidth',2 );
        plot((1:max_iter), (est_data.h_error),'-r*','linewidth',2 );
        %title('\fontsize{30}L^2 error in h');
        xlabel('\fontsize{30}number of iterations');
        ylabel('\fontsize{30}L2 error');
        ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
        %str = sprintf('L^2 error: %0.2e',est_data.h_error(max_iter));
        %temp_text = text(max_iter/3, 0.5*min(log10(est_data.h_error(2:end))),str);
        %set(temp_text,'fontsize',32);
                
    elseif true_data.use_data == 0
        fprintf('est. sum square L2norm %e\n',est_data.final_L2cost);
        fprintf('est ortho. %e\n', est_data.ortho(max_iter));
    end
      
    fig_mean_f = figure(12);
    plot(t,true_data.mean_f,'linewidth',2); My_Figure(7,6);
    xlim([t(1) t(end)]);
    xlabel('\fontsize{30}time');
    ylabel('\fontsize{30}value');
    %title('\fontsize{30}mean of observed signals f_{i}(t)');
    ax = gca; ax.XAxis.FontSize = 20; ax.YAxis.FontSize = 20;
    fprintf('time cost is %f seconds\n',est_data.timecost);
    
%%%%%  save results: publish and save data and pictures  %%%%%%        
    temp_dir = pwd;
    if( isdir('html') )
        rmdir('html','s');
    end
    mkdir('html');
    cd('html');
    
    print(fig_f,'-dpng','f');    
    if true_data.use_data == 1 
        print(fig_g,'-dpng','g');
        print(fig_h,'-dpng','h');
        print(fig_r,'-dpng','r');
        
        print(fig_L2cost_est_g,'-dpng','L2cost_est_g');
        print(fig_L2cost_est_h,'-dpng','L2cost_est_h');
        print(fig_L2cost_gam,'-dpng','L2cost_gam');
    end    
    
    print(fig_est_g,'-dpng','est_g');
    print(fig_est_h,'-dpng','est_h');
    print(fig_gam,'-dpng','gam');
    print(fig_L2cost,'-dpng','L2cost');

    print(fig_mean_f,'-dpng','mean_f');  
    save('F_data.mat','para','true_data','est_data');
    
    %%%%%%%%
    fid = fopen('Result.txt','wt');
    fprintf(fid,'implement time:%d.%d.%d.%d.%d.%d\n\n',fix(clock));
    fprintf(fid, 'T: %d, N: %d\n',true_data.T,true_data.N);
    
    if true_data.use_data == 1 % synthetic data
        fprintf(fid,'\ntrue sum square L2norm %e\n',est_data.syn_L2cost);
        fprintf(fid,'est. sum square L2norm %e\n',est_data.final_L2cost);
        fprintf(fid,'abslote L2_error of r %e\n', est_data.r_error(max_iter));
        fprintf(fid,'relative L2_error of g %e\n', est_data.g_error(max_iter) );
        fprintf(fid,'relative L2_error of h %e\n\n', est_data.h_error(max_iter) );
        fprintf(fid,'true ortho. %e\n', est_data.true_ortho );
        fprintf(fid,'est ortho. %e\n', est_data.ortho(max_iter)); 
    end
    fclose(fid);
    
    cd(temp_dir); close all;  
end