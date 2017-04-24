% this is the main function for trend separation, 
% continued on old version which is created on July, 27, 2014
% created on March, 21, 2015
% est_energy for orthogonal penalty to be ajusted, July, 23, 2015
% may use Simpson's rule to replace evaluating a quadrature!!
%            current version is using spline, may to too expansive
%                                            Aug. 12, 2015
% continued on Sep. 7, 2015, add change directory and save data, pictures
% continued on July, 2, 2016; this is jusing ploting results on the screen                           

function Plot_Result_on_Screen(para, true_data, est_data)

    max_iter = para.max_iter;
    t = true_data.t;

    if true_data.use_data == 1 %% synthetic       
        fig1 = figure(1);
        subplot(2,2,1); 
        plot(t,true_data.g);
        title('true g'); xlabel('time'); ylabel('value');
    
        subplot(2,2,2);
        plot(t,true_data.r);
        title('random warping functions');
        
        subplot(2,2,3);
        plot(t,true_data.h);
        title('true trend h'); xlabel('time'); ylabel('value');
        
        subplot(2,2,4);
        plot(t,true_data.f);
        title('observed signal f_{i}(t)'); xlabel('time'); ylabel('value');  
        
    elseif true_data.use_data == 0 %% real or import data
        fig1 = figure(1);
        plot(t,true_data.f);
        title('observed signal f_{i}(t)'); xlabel('time'); ylabel('value');
    end
    
%%%%%%%%%%%   plot estimated solutions   %%%%%%%%%%%%%%%    
    est_g = est_data.est_g;
    est_h = est_data.est_h;
    gam = est_data.gam;
   
    fig2 = figure(2); 
    subplot(2,2,1);  
    if true_data.use_data == 0
        plot(t,est_g(:,max_iter),'--r');
        legend('est g','location','Best');
        title('est g');
        
    elseif true_data.use_data == 1
        plot(t,true_data.g,t,est_g(:,max_iter),'--r');
        legend('true g','est g','location','Best');
        title('true g and est g');
    end
    xlabel('time'); ylabel('value');
   
    subplot(2,2,2);
    plot(t,gam(:,:,max_iter));
    title('est. gamma');
    
    subplot(2,2,3);
    if true_data.use_data == 0
        plot(t,est_h(:,max_iter),'--r'); 
        legend('est h','location','Best');
        title('est h');
        
    elseif true_data.use_data == 1
        plot(t,true_data.h,t,est_h(:,max_iter),'--r');
        legend('true h','est h','location','Best');
        title('true h and est h');
    end
    xlabel('time'); ylabel('value');

%%%%%%%%%   plot error   %%%%%%%%%%
    subplot(2,2,4);
    plot((1:max_iter), log10(est_data.L2cost(3,:)),'r'); hold on;
    plot((1:max_iter), log10(est_data.L2cost(3,:)),'r*');
    xlabel('# of iterations');
    ylabel('log10 square L^2');
    title('iteration of log10 square L^2');

    if true_data.use_data == 1
        fprintf('est. sum square L2norm %e\n',est_data.L2cost(3,max_iter));
        fprintf('abslote L2_error of r %e\n', est_data.r_error(max_iter));
        fprintf('relative L2_error of g %e\n', est_data.g_error(max_iter));
        fprintf('relative L2_error of h %e\n\n', est_data.h_error(max_iter));
        fprintf('true ortho. %e\n', est_data.true_ortho );
        fprintf('est ortho. %e\n', est_data.ortho(max_iter));

        fig3 = figure(3);
        subplot(2,2,1);
        plot((1:max_iter), log10(est_data.r_error) );
        title('Log10 L^2 error in r_i');
        xlabel('iter');

        subplot(2,2,2);
        plot((1:max_iter), log10(est_data.g_error) );
        title('Log10 L^2 error in g');
        xlabel('iter');

        subplot(2,2,3);
        plot((1:max_iter), log10(est_data.h_error) );
        title('Log10 L^2 error in h');
        xlabel('iter');
        
    elseif true_data.use_data == 0
        fprintf('est. sum square L2norm %e\n',est_data.final_L2cost); 
        fprintf('est ortho. %e\n', est_data.ortho(max_iter));
    end
    
    fig4 = figure(4);
    plot(t,true_data.f);
    xlabel('time');
    ylabel('value');
    title('observed signals f_{i}(t)');
    
    fig5 = figure(5);
    plot(t, true_data.mean_f);
    xlabel('time'); ylabel('value'); title('mean of observed signals f_{i}(t)');
   
    fprintf('time cost is %f seconds\n',est_data.timecost);
        
%%%%%  save results: publish and save data and pictures  %%%%%%        
    root_dir = pwd;
    if( isdir('html') )
        rmdir('html','s');
    end
    mkdir('html');
    cd('html');
    print(fig1,'-dpng','a');
    print(fig2,'-dpng','b');
    if true_data.use_data == 1, print(fig3,'-dpng','c'); end
    print(fig4,'-dpng','observed_f');
    print(fig5,'-dpng','mean_f');
    
    save('F_data.mat','para','true_data','est_data');
    cd(root_dir);   
end