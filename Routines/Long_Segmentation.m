% this is a program for long signal segmentation
% idea from Kurtek 2013 paper.
% created at April, 9 , 2015
% continued on April 23, 2015

clear; close all;
my_import;

display = 0;  % 0 is no display; 1 is to display 

% for all ecg data
load('ecgca102_edfm.mat','val');
f_origin = val';
[T_origin,N_origin] = size(f_origin);
T_grid = (1:T_origin)';

% plot original results
if display == 1
    figure(1);
    for i = 1:N_origin-2
        subplot(2,2,i);
        plot((1:T_origin),f_origin(:,i));
    end
end

% downsize by 1/10
temp_grid = floor(linspace(1,T_origin,T_origin/10))';
f = zeros(length(temp_grid),N_origin);

for i=1:N_origin-2
    f(:,i) = spline(T_grid,f_origin(:,i),temp_grid);
end

if display == 1
    figure(2);
    for i = 1:N_origin-2
        subplot(2,2,i)
        plot(temp_grid,f(:,i));
    end
    
    figure(3);
    plot(T_grid,f_origin(:,1),temp_grid,f(:,1));
    legend('original','splined');
    title('difference between original and simplified');
end

%%%%%%%  set template manually  %%%%%%%
t_grid = (1:length(temp_grid));
figure(4);
plot(t_grid,f(:,1));

win_start = 45;
win_end = 125;
win_length = win_end - win_start;

increment = 10;
iterations = floor( (length(t_grid)-win_length) /increment)-1;

f_template = f((win_start:win_end),1);
t = linspace(0,1,win_length+1);
q_f_template = SRVF(t,f_template);

gam = zeros(win_length+1,iterations);     % initialized with indentity


for i = 1:iterations
    f_window = f( 1 + i*increment : win_length + 1 + i*increment,1);
    q_f_window = SRVF(t,f_window);
    
    % do DynamicProgramming
    lambda = 0;
    gam_temp = DynamicProgrammingQ(q_f_window',q_f_template',lambda,0);
    gam(:,i) = (gam_temp-gam_temp(1))/(gam_temp(end)-gam_temp(1));
    
    fprintf('finish %d iterations\n',i);
    
    % do transformation to get phi
    
end
