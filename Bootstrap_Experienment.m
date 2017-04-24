%%%%%%%%%% generate bootstrap and hypothesis testing, confidence interval
% created at June, 29,2016; 
% continued on July,3,2016; finished for generating bootstrap :)
% continued on July,6,2016; doing hypothesis testing
% continued on July,11,2016; output results to filess
% continued on July,16,2016; all routine for h now can for g! done;
%%  input generate bootstrap replicates  %%%%%%%%%%%

clear; close all; 
fprintf('implement time:%d.%d.%d.%d.%d.%d\n',fix(clock));  % display the time
my_import;

Bootstrap_Number = 1;
max_iter = 10;         % in each iteration

% load data
[syn_para] = Set_Synthetic_Data_Para;              fprintf(' finish Set_Synthetic_Data_Para \n');
[true_data] = Generate_Synthetic_Data(syn_para);   fprintf(' finish Generate_Synthetic_Data \n');
%[true_data] = Import_Data;                        fprintf('\n finish Import_Data \n');

est_g = zeros(true_data.T,Bootstrap_Number);
est_h = zeros(true_data.T,Bootstrap_Number);
gam = zeros(true_data.T,true_data.N,Bootstrap_Number);
L2cost = zeros(max_iter,Bootstrap_Number);
final_L2cost = zeros(1,Bootstrap_Number);
ortho = zeros(1,Bootstrap_Number);
    
for B = 1:Bootstrap_Number
fprintf('doing %d bootstrap... \n',B);
% random selection with replacement
 Bootstrap = datasample(true_data.f,true_data.N,2);
 
% run Trend extraction with all basis
    min_h_index = 1; % if start from 2, then add 1 when input to Set_Algorithm_Para
    max_h_index = 10;
    
    for basis_type = 1:4                % from 1 to 4 are Cos, Sin, MFS,sLegen
        for basis_number = min_h_index : max_h_index  % from 1,2 to 10

            [bootstrap_true_data] = Import_Data(Bootstrap);          
            fprintf('\n finish Import_Data \n');
            [para] = Set_Algorithm_Para(max_iter, 10+5*basis_type ,min_h_index,basis_number);  
            fprintf(' finish Set_Algorithm_Para \n');
            [est_data2(basis_number,basis_type)] = Trend_Separation(para, bootstrap_true_data);
            fprintf('\n finish Trend_separation \n');
            
        end %%%%%%%%%% end of multiple iteration of h index
    end %%%%%%%%%% end of multiple h basis
      
    % pick up the basis type and index that has lowest L^2 cost
    temp2 = cell2mat({est_data2.final_L2cost}); 
    [~,idx2] = min(temp2);
    idx2_basis = ceil(idx2/(max_h_index-min_h_index+1));
    idx2_number = idx2 - (max_h_index-min_h_index+1)*(idx2_basis - 1);
     
    % pick up desired results
    est_g(:,B) = est_data2(idx2_number,idx2_basis).est_g(:,max_iter);
    est_h(:,B) = est_data2(idx2_number,idx2_basis).est_h(:,max_iter); 
    gam(:,:,B) = est_data2(idx2_number,idx2_basis).gam(:,:,max_iter);
    L2cost(:,B) = est_data2(idx2_number,idx2_basis).L2cost(3,:)';
    final_L2cost(B) = est_data2(idx2_number,idx2_basis).L2cost(3,max_iter);
    ortho(B) = est_data2(idx2_number,idx2_basis).ortho(:,max_iter);
    
end % end of Bootstrap loops

    save('bootstrap_replicate.mat','est_g','est_h','gam','L2cost','final_L2cost','ortho','true_data');   
    
%%  run hypothesis testing and confidence interval  %%%%%%%%%%%
    clear;
    %load bootstrap_replicate.mat;
    %load bootstrap_replicate_syn_500_N10;
    load bootstrap_replicate_electricity_500;
    
    mu = 0; sigma = 1; p = 0.975; z_score = norminv(p,mu,sigma);
    
    t = true_data.t; [T,B] = size(est_h);
    replicate_h_0 = zeros(1,B);
    replicate_h_c = zeros(1,B);
    replicate_g_0 = zeros(1,B);
    replicate_g_c = zeros(1,B);
    
    for b = 1:B
        replicate_h_0(b) = L2norm(t, est_h(:,b));
        avg_h = Simpsons_Rule(est_h(:,b),0,1);
        replicate_h_c(b) = L2norm(t,est_h(:,b)-avg_h);
        
        replicate_g_0(b) = L2norm(t, est_g(:,b));
        avg_g = Simpsons_Rule(est_g(:,b),0,1);
        replicate_g_c(b) = L2norm(t,est_g(:,b)-avg_g);
    end
    
    % testing zero and constant for h and g
    % testing zero and constant for h and g  
    se_replicate_h_0 = Bootstrap_se(replicate_h_0);
    theta_h_0 = L2norm(t,true_est_h);
    
    se_replicate_h_c = Bootstrap_se(replicate_h_c);
    theta_h_c = L2norm(t,true_est_h - Simpsons_Rule(true_est_h,0,1));
    
    se_replicate_g_0 = Bootstrap_se(replicate_g_0);
    theta_g_0 = L2norm(t,true_est_g);
    
    se_replicate_g_c = Bootstrap_se(replicate_g_c);
    theta_g_c = L2norm(t,true_est_g - Simpsons_Rule(true_est_g,0,1));
    
    fprintf('H0:zero est h; p-value: %e; ',1-cdf('Normal',theta_h_0/se_replicate_h_0,mu,sigma));
    fprintf('C.I. %f to %f\n',theta_h_0-z_score*se_replicate_h_0,theta_h_0+z_score*se_replicate_h_0);

    fprintf('H0:constant est h; p-value: %e; ',1-cdf('Normal',1-theta_h_c/se_replicate_h_c,mu,sigma));
    fprintf('C.I. %f to %f\n', theta_h_c - z_score*se_replicate_h_c, theta_h_c + z_score*se_replicate_h_c);
    
    fprintf('H0:zero est g; p-value: %e; ',1-cdf('Normal',theta_g_0/se_replicate_g_0,mu,sigma));
    fprintf('C.I. %f to %f\n',theta_g_0-z_score*se_replicate_g_0,theta_g_0+z_score*se_replicate_g_0);

    fprintf('H0:constant est g; p-value: %e; ',1-cdf('Normal',1-theta_g_c/se_replicate_g_c,mu,sigma));
    fprintf('C.I. %f to %f\n', theta_g_c - z_score*se_replicate_g_c, theta_g_c + z_score*se_replicate_g_c);
    
    % pointwise confidence interval for h and g
    [ upper_est_h,lower_est_h ] = Confidence_Interval( true_est_h, est_h );
    [ upper_est_g,lower_est_g ] = Confidence_Interval( true_est_g, est_g );

    %%%%%%%%%%%%%%%%%%%%%%%%% plotting here %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% plotting here %%%%%%%%%%%%
    fig_replicate_h = figure(1); My_Figure(7,6);
    plot(t,est_h,'linewidth',2);
    ylim([-3 3]);
    xlabel('time'); ylabel('value');
    title('Replicates of est. h');
    
    fig_replicate_g = figure(2); My_Figure(7,6);
    plot(t,est_g,'linewidth',2);
    ylim([-3 3]);
    xlabel('time'); ylabel('value');
    title('Replicates of est. g');
    
    fig_replicate_h_0 = figure(3); My_Figure(7,6);
    hist(replicate_h_0); 
    xlabel('L2norm of h(t)'); ylabel('frequency');
    title('histogram of L2norm of h(t)');
    
    fig_replicate_h_c = figure(4); My_Figure(7,6);
    hist(replicate_h_c); 
    xlabel('L2norm of (h-avg(h))'); ylabel('frequency');
    title('histogram of L2norm of (h-avg(h))');
    
    fig_replicate_g_0 = figure(5); My_Figure(7,6);
    hist(replicate_g_0); 
    xlabel('L2norm of g(t)'); ylabel('frequency');
    title('histogram of L2norm of g(t)');
    
    fig_replicate_g_c = figure(6); My_Figure(7,6);
    hist(replicate_g_c); 
    xlabel('L2norm of (g-avg(g))'); ylabel('frequency');
    title('histogram of L2norm of (g-avg(g))');
    
    fig_conf_band_h = figure(7); My_Figure(7,6);
    plot(t,true_est_h,t,upper_est_h,t,lower_est_h,'linewidth',2);
    ylim([-3 3]);
    xlabel('time'); ylabel('value');
    title('confidence band for the trend h(t)');
    legend('est. h','upper bound','lower bound','location','northeastoutside');
    
    fig_conf_band_g = figure(8); My_Figure(7,6);
    plot(t,true_est_g,t,upper_est_g,t,lower_est_g,'linewidth',2);
    ylim([-3 3]);
    xlabel('time'); ylabel('value');
    title('confidence band for main structure g(t)');
    legend('est. g','upper bound','lower bound','location','best');
    
    % save figures % save figures % save figures % save figures % save figures
    root_dir = pwd;
    if ( isdir('bs_plots') );
        rmdir('bs_plots','s');
    end  
    mkdir('bs_plots'); cd('bs_plots');
    
    print(fig_replicate_h,'-dpng','replicate_h');
    print(fig_replicate_g,'-dpng','replicate_g');
    print(fig_replicate_h_0,'-dpng','replicate_h_0');
    print(fig_replicate_h_c,'-dpng','replicate_h_c');
    print(fig_replicate_g_0,'-dpng','replicate_g_0');
    print(fig_replicate_g_c,'-dpng','replicate_g_c');
    print(fig_conf_band_h,'-dpng','conf_band_h');
    print(fig_conf_band_g,'-dpng','conf_band_g');
    cd(root_dir); close all;