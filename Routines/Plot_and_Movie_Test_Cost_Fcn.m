% this routine plots "Test_Cost_Fcn.m"
% creat on Nov. 2, 2014; successful

clear; close all;

load F_data.mat;

t = true_data.t;
%g = true_data.g;
%h = true_data.h;
f = true_data.f;
%r = true_data.r;
%energy = error.energy;

gam = est_data.gam;
est_g = est_data.est_g;
est_h = est_data.est_h;
est_energy = error.est_energy;


display = 0;  % 1 is plot; 0 is no plot
movie = 1;    % movie 

%stop = j;
stop = 50;
energy_plot = 2;

%%%%%%%%%%%%%% start %%%%%%%%%%%%%%
if display == 1;
    figure(1);
    subplot(2,2,1);
    plot(t,g);
    title('true g(t)');
    xlabel('time');
    ylabel('value');
    
    subplot(2,2,2);
    plot(t,r);
    title('random warping functions');
    
    subplot(2,2,3);
    plot(t,h)
    title('true trend h');
    xlabel('time');
    ylabel('value');
    
    subplot(2,2,4);
    plot(t,f);
    title('observed signal \{f_i\}');
    xlabel('time');
    ylabel('value');
    
    
    %%%%%%%%%%%%%% middle %%%%%%%%%%%%%%
    
    for j = 1:stop
        
        figure(2);
        subplot(2,2,2);
        plot(t,gam(:,:,j)); title('est. gamma');
        
        subplot(2,2,3);
        plot(t,h,t,est_h(:,j),'--r');
        xlabel('time'); ylabel('value');
        legend('true h','est h','location','Best');
        title('true h and est h');
        
        subplot(2,2,1);
        plot(t,g,t,est_g(:,j),'--r');
        xlabel('time'); ylabel('value');
        legend('true g','est g','location','Best');
        title('true g and est g');
        
        subplot(2,2,4);
        tempT = (1:stop);
        if energy_plot == 1
            plot(tempT,log(energy(1:stop)),'b',tempT,log(est_energy(1:stop)),'r:');
            legend('true energy','est. energy','location','Best');
            hold on;
            
            plot(tempT,energy(1:stop),'bo',tempT,est_energy(1:stop),'r*');
            xlabel('# of iterations');
            ylabel('minimized energy');
            title('square L^2: true energy and est. energy');
            
        elseif energy_plot == 2
            % plot(tempT,log(est_energy(1:stop)),'r'); hold on;
            % plot(tempT,log(est_energy(1:stop)),'r*');
            
            plot(j,log(est_energy(j)),'r'); hold on;
            plot(j,log(est_energy(j)),'r*');
            
            legend('est. energy','location','Best');
            xlabel('# of iterations');
            ylabel('minimized energy');
            title('log square L^2 of est. energy');
        end
    end
    
    %%%%%%%%% end %%%%%%%%%%%%%%
    
    
    figure(3);
    subplot(1,3,1);
    plot(tempT,gam_error(1:stop));
    title('L^2 error in r_i');
    xlabel('iter');
    
    subplot(1,3,2);
    plot(tempT,g_error(1:stop));
    title('L^2 error in g');
    xlabel('iter');
    
    subplot(1,3,3);
    plot(tempT,h_error(1:stop));
    title('L^2 error in h');
    xlabel('iter');
end

%%%%%%%%%%%%%%%%%%% movie %%%%%%%%%%%%%%%%%%%5
if movie == 1
    writerObj = VideoWriter('sine.avi');
    writerObj.FrameRate = 5; % vedio speed ! % big number is faster
    open(writerObj);
 
    
    for j = 1:stop
        
        figure(2);
        
        subplot(2,3,1);
        plot(t,f);
        title('observed signal \{f_i\}');
        xlabel('time');
        ylabel('value');
        
        subplot(2,3,2);
        plot(t,g,t,est_g(:,j),'--r');
        xlabel('time'); ylabel('value');
        legend('true g','est g','location','Best');
        title('true g and est g');
        
        subplot(2,3,3);
        plot(t,gam(:,:,j)); title('est. gamma');
        
        subplot(2,3,5);
        plot(t,h,t,est_h(:,j),'--r');
        xlabel('time'); ylabel('value');
        legend('true h','est h','location','Best');
        title('true h and est h');
        
        
        subplot(2,3,6);
        tempT = (1:stop);
        if energy_plot == 1
            plot(tempT,log(energy(1:stop)),'b',tempT,log(est_energy(1:stop)),'r:');
            legend('true energy','est. energy','location','Best');
            hold on;
            
            plot(tempT,energy(1:stop),'bo',tempT,est_energy(1:stop),'r*');
            xlabel('# of iterations');
            ylabel('minimized energy');
            title('square L^2: true energy and est. energy');
            
        elseif energy_plot == 2
            % plot(tempT,log(est_energy(1:stop)),'r'); hold on;
            % plot(tempT,log(est_energy(1:stop)),'r*');
            
            plot(j,log(est_energy(j)),'r'); hold on;
            plot(j,log(est_energy(j)),'r*');
            
            legend('est. energy','location','Best');
            xlabel('# of iterations');
            ylabel('minimized energy');
            title('log square L^2 of est. energy');
        end
        
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        
        
    end
    
    close(writerObj);
end