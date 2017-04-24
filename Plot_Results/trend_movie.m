% this is a program for producing movie
% created on Oct. 14, 2016
clear; close all; load F_data;

 writerObj = VideoWriter('movie.avi');
 writerObj.FrameRate = 1;
 open(writerObj);

fig=figure(1); %set(gcf, 'Color', default);
set(fig, 'Position', [-120 100 1300 300])
for iter = 1:10
    subplot(1,4,1); 
    plot(true_data.t,true_data.f,'linewidth',2); 
    title('observed signal f_{i}(t)','FontSize', 19); 
    xlabel('time','FontSize', 14); ylabel('value','FontSize', 14);  
    xlim([0,1]);
        
    subplot(1,4,2);
    plot(true_data.t,true_data.h,'b', true_data.t,est_data.est_h(:,iter),'--r','linewidth',2); 
    legend('true h','est h','location','northeast');
    xlabel('time','FontSize', 14); ylabel('value','FontSize', 14);  
    title('Trend Estimation','FontSize', 19);
    ylim([-0.5,1.5]); xlim([0,1]);
    
    subplot(1,4,3);
    plot(true_data.t,true_data.g,'b',true_data.t,est_data.est_g(:,iter),'--r','linewidth',2); 
    legend('true g','est g','location','northeast');
    xlabel('time','FontSize', 14); ylabel('value','FontSize', 14);  
    title('Seasonality Estimation','FontSize', 18);
    ylim([-1.5,1.5]); xlim([0,1]);
     
    subplot(1,4,4);
    plot((1:10),(-4:4/9:0),'w', (1:iter),log10(est_data.Check_Cost(4,(1:iter))),'r*'); 
    xlabel('# of iterations','FontSize', 14);
    ylabel('log10 L^2 cost','FontSize', 13);
    title('Minimization of cost function','FontSize', 15);
    xlim([0,10]);
    
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);

%%
% writerObj = VideoWriter('sine.avi');
% writerObj.FrameRate = 10;
% open(writerObj);
% 
% n = 100;
% t = linspace(0,1,n);
% y = sin(2*pi*t);
% 
% plot(t,y,'linewidth',2); hold on;
% for i=1:n
%     plot(t(i),y(i),'r*','linewidth',5);
%     frame = getframe;
%     writeVideo(writerObj,frame);
% end
% 
% close(writerObj);

%%%% tips for creating movie
% open object
% for j = 1:n 
%    plotting commands
%    F(j) = getframe;
% end
% movie(F)
% close object
