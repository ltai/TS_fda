%%%%%%%% this is an old version of comparing all models.
% However, only one iteration. this file can be deleted
% on Feb. 4, 2016

% this is a routine for comparing other simple trend spearation model
% all call existing routines
%
% created on Sep. 23, 2015; successful!
%                        let outputing name adjusted
% continued on sep. 28, 2015; 
% continued on Oct. 12, 2015; add all relative error of g,h and gam, 
% continued on Nov. 4, 2015; all model change to multiple iterations

clear; 
close all;
fprintf('implement time:%d.%d.%d.%d.%d.%d\n ',fix(clock));  % display the time
%my_import; 
tic;

model_choice = 2;   % 1 is cross-sectional;                               
                    % 2 is projection (no warping);          
                    % 3 is alignment only;  
                    % 4 is projection then alignment;
                    % 5 is alignment then projection;    
                    % 6 is signal separation;

my_data_choice = 1;    % 1 is artificial, simulated here; 
                       % 2 is real data  
                       % 3 is artifical, load from F_data.

display_figure = 2;  % 1 is not show; 2 is show
        

%%%%% data input %%%%%
if my_data_choice == 1             %%%%% artificial data
    T = 101; N = 5;
    
    t = linspace(0,1,T)';
    fcn = inline('sin(2*pi*t)','t'); % 5pi is same as paper
    true_data.g = fcn(t);
    true_data.h = 0*t + 1;
    true_data.r = get_gam(T,N,4);
    true_data.r = Gamma_with_inverse_kmean_id(true_data.r);
    % r = Gamma_with_inverse_Amean_id(r);
    true_data.f = zeros(T,N);
    
    for i = 1:N
        true_data.f(:,i) = MyGroupAction( t, true_data.g, true_data.r(:,i), 3) + exp(-t);
    end
    f = true_data.f;
    
    para.h_start_idx = 1;
    para.h_end_idx = 4;
    max_iter = 20;
    
elseif my_data_choice == 2         %%%% real data: growth velocity
     load F_data; 
     
     f = true_data.f;
     [T,N] = size(f);
     t = linspace(0,1,T)';
    
elseif my_data_choice == 3  %%% artificial data, but load from F_data
    load F_data;  
    
    f = true_data.f;
    [T,N] = size(f);
    t = linspace(0,1,T)';
end

%%%%% choose method %%%%%
%%%%% choose method %%%%%
if(true)
    if model_choice == 1             % 1 is cross-sectional; 
        
        est_g = zeros(T,max_iter);
        L2cost = zeros(max_iter,1);
        
        for j = 1:max_iter
            est_g(:,j) = mean(f,2);
        end
        
        for j = 1:max_iter
            temp_L2cost = 0;
            for i = 1:N
                temp_L2cost = temp_L2cost + SquareL2norm(t,f(:,i) - est_g(:,j));
            end
            L2cost(j) = temp_L2cost/N;
        end
        
        all_est_g = est_g; 
        est_g = est_g(:,end);
        
    elseif model_choice == 2         % 2 is projection (no warping), para. h!
        h_start_idx = 1;
        h_end_idx = 4;
        
        max_iter = 20;
        est_g = zeros(T,max_iter);
        est_h = zeros(T,max_iter);
        L2cost = zeros(max_iter,1);

        est_g(:,1) = mean(f,2);
        
        for j = 2:max_iter
            % update est g
            temp_vec = zeros(T,1);
            for i = 1:N
                temp_vec = temp_vec + f(:,i) - est_h(:,j-1);
            end
            est_g(:,j) = temp_vec/N;
            
            proj_g = L2_Shift_Legendre_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
            est_g(:,j) = est_g(:,j) - proj_g;
            
            % update est h
            temp_vec = zeros(T,1);
            for i = 1:N
                temp_vec = temp_vec + f(:,i) - est_g(:,j);
            end
            est_h(:,j) = temp_vec/N;   
            est_h(:,j) = L2_Shift_Legendre_Approximation(est_h(:,j), h_start_idx, h_end_idx);

        end
       
        for j = 1:max_iter
            temp_L2cost = 0;
            for i = 1:N
                temp_L2cost = temp_L2cost + SquareL2norm(t,f(:,i) - est_g(:,j) - est_h(:,j) );
            end
            L2cost(j) = temp_L2cost/N;
        end
        
        all_est_g = est_g; est_g = est_g(:,end);
        all_est_h = est_h; est_h = est_h(:,end);

    elseif model_choice == 3  % 5 is alignment only
        [fn,gam,qn,mq,gamI] = mainWarpingWrapper(t,f,0);
        est_g = mean(fn,2);
        
        L2cost = 0;
        for i = 1:N
            L2cost = L2cost + SquareL2norm(t,f(:,i) - MyGroupAction(t,est_g,gam(:,i),1));
        end
        L2cost = L2cost/N;
        
    elseif model_choice == 4         % 4 is projection then alignment 
        %%% projection to get h
        h_start_idx = 1;
        h_end_idx = 4;
        temp_g = mean(f,2);
        est_h = L2_Shift_Legendre_Approximation(temp_g, h_start_idx, h_end_idx);

        %%% subtraction to get f_i - h in order to get g_i
        temp_f = zeros(T,N);
        for i = 1:N
            temp_f(:,i) = f(:,i) - est_h;    
        end

        %%% alignment to get g
        [fn,gam,qn,mq,gamI] = mainWarpingWrapper(t,temp_f,0);   
        est_g = mean(fn,2);

        L2cost = 0;
        for i = 1:N
            L2cost = L2cost + SquareL2norm(t,f(:,i) - MyGroupAction(t,est_g,gam(:,i),1) - est_h);
        end
        L2cost = L2cost/N;

    elseif model_choice == 5         % 5 is alignment then projection
        %%% alignment
        [fn,gam,qn,mq,gamI] = mainWarpingWrapper(t,f,0);
        temp_g = mean(fn,2);

        %%% projection
        h_start_idx = 1;
        h_end_idx = 4;

        est_h = L2_Shift_Legendre_Approximation(temp_g, h_start_idx, h_end_idx);
        est_g = temp_g - est_h;

        L2cost = 0;
        for i = 1:N
            L2cost = L2cost + SquareL2norm(t,f(:,i) - MyGroupAction(t,est_g + est_h,gam(:,i),1));
        end
        L2cost = L2cost/N;

    end
end

%%%%%%%%%  compute energy, for artificial data only %%%%%%%%%%%%%
if my_data_choice == 3 || my_data_choice == 1
   if(exist('est_g','var'))
       g_error = L2norm(t, est_g - true_data.g);
       fprintf('\nrelative L2_error of g %e\n', g_error/L2norm(t,true_data.g) );
   end
   
   if(exist('est_h','var')) 
       h_error = L2norm(t, est_h - true_data.h); 
       fprintf('relative L2_error of h %e\n', h_error/L2norm(t,true_data.h) );
   end
   
   if(exist('gam','var'))
       gam_error = 0;
       for i = 1:N
           gam_error = gam_error + L2norm(t, gam(:,i) - true_data.r(:,i));
       end
       gam_error = gam_error/N;
       fprintf('abslote L2_error of r %e\n', gam_error);
   end
end

%%%%%%%%%%%%% check orthogonality %%%%%%%%%%%
if exist('est_g','var') && exist('est_h','var')
    if my_data_choice == 3 || my_data_choice == 1
        fprintf('true ortho. %e\n', Check_Orthogonality(true_data.g, true_data.h, t) );
        fprintf('est ortho. %e\n', Check_Orthogonality(est_g(:,end), est_h(:,end), t));
    end
end

%%%%% plot and save %%%%%
%%%%% plot and save %%%%%
if display_figure == 1 
    set(0,'DefaultFigureVisible','off');
elseif display_figure == 2 
    set(0,'DefaultFigureVisible','on');
end
warning('off', 'all');   %%% surpress warning

%%%%%  plot f
fig_f = figure(1);  My_Figure(7,6); 
if my_data_choice == 3 || my_data_choice == 1
    plot(t,f,'linewidth',2); 
    title('\fontsize{36}observed signals {f_i}(t)');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 
    
elseif my_data_choice == 2
    plot(linspace(0,18,T)',f,'linewidth',2); 
    xlim([0 18]); 
    title('\fontsize{36}male growth velocity f_{i}(t)');
    xlabel('\fontsize{36}age(years)');
    ylabel('\fontsize{30}growth velocity(cm/year)');
end

%%%%% plot est g
fig_est_g = figure(2); My_Figure(7,6);

if my_data_choice == 3 || my_data_choice == 1 
    plot(t,est_g,'--r','linewidth',2);
    hold on;
    plot(t,true_data.g,'b','linewidth',2);
    legend('\fontsize{26}est. g','\fontsize{26}g','location','best');
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value');
    
elseif my_data_choice == 2
    plot(linspace(0,18,T)',est_g,'--r','linewidth',2);
    xlim([0 18]); 
    title('\fontsize{36}est. g');
    xlabel('\fontsize{36}age(years)');
    ylabel('\fontsize{30}growth velocity(cm/year)');
    
end

%%%%%% plot est h
if model_choice == 2 || model_choice == 4 || model_choice == 5 ...
                                          || model_choice == 6 %% have h
    fig_est_h = figure(3); My_Figure(7,6);
    if my_data_choice == 3 || my_data_choice == 1
        plot(t,est_h,'--r','linewidth',2);
        hold on;
        plot(t,true_data.h,'b','linewidth',2);
        legend('\fontsize{26}est. h','\fontsize{26}h','location','best');
        title('\fontsize{36}est. h');
        xlabel('\fontsize{36}time');
        ylabel('\fontsize{36}value'); 
        
    elseif my_data_choice == 2
        plot(linspace(0,18,T)',est_h,'--r','linewidth',2);
        xlim([0 18]); 
        title('\fontsize{36}est. h');
        xlabel('\fontsize{36}age(years)');
        ylabel('\fontsize{30}growth velocity(cm/year)');
    end

end   

%%%% plot gamma %%%%
if model_choice == 3 || model_choice == 4 || model_choice == 5 ...
                                          || model_choice == 6 %% have gamma
    fig_gam = figure(4); My_Figure(5.8,6);
    plot(linspace(0,1,T)',gam,'linewidth',2);
    title('\fontsize{36}est. gamma');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}time'); 
end

if model_choice == 3
    fig_aligned_f = figure(5); My_Figure(7,6);
    %plot(t,fn,'linewidth',2);
    %title('\fontsize{36}alignmed signals f_i(t)');
    %xlabel('\fontsize{36}time');
    %ylabel('\fontsize{36}value'); 
    
    plot(linspace(0,18,T)',fn,'linewidth',2);
    xlim([0 18]); 
    title('\fontsize{36}alignmed growth data f_i(t)');
    xlabel('\fontsize{36}age(years)');
    ylabel('\fontsize{30}growth velocity(cm/year)');
end

if model_choice == 1
    fig_mean_f = figure(6); My_Figure(7,6);
    plot(t,mean(f,2),'linewidth',2);
    title('\fontsize{36}mean of observed f_i(t)');
    xlabel('\fontsize{36}time');
    ylabel('\fontsize{36}value'); 
end

fprintf('est. sum square L2norm %e\n',L2cost);
        
toc;

%%%%%%%%%%%%%%%  output data, save and print figure %%%%%%%%%%%%%%%%%%%%%
    root_dir = pwd;
    if( isdir('html') )
        rmdir('html','s');
    end
    mkdir('html');
    cd('html');
    print(fig_f,'-dpng','f');
    print(fig_est_g,'-dpng','est_g');
    
    if model_choice == 2 || model_choice == 4 || model_choice == 5 ...
                                          || model_choice == 6 %% have h
        print(fig_est_h,'-dpng','est_h');
    end    
    if model_choice == 3 || model_choice == 4 || model_choice == 5 ...
                                          || model_choice == 6 %% have gamma
        print(fig_gam,'-dpng','gam');
    end
    if model_choice == 3
        print(fig_aligned_f,'-dpng','aligned_f');
    end
    if model_choice == 1
        print(fig_mean_f,'-dpng','mean_f');
    end
    save('F_data.mat');
    cd(root_dir);   
    %close all;