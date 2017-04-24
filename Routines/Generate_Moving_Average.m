% this is a test for moving average method.
% created on September 14, 2015

clear;
close all;

method_idx = 2; % 1 is filter; 2 is convolution
windowWidth = 3; 

% data = (1:20)'+ 5*randn(20,1);
% load TLH_AverageTemperature; data = f;
% load Florida_Sea_Level; data = f;
 load US_Electricity_Price_Difference_2005_2010; data = f;

[T,N] = size(data);
grid = (1:T)';
output_f = zeros(T,N);

if method_idx == 1
    kernel = ones(windowWidth,1) / windowWidth;
    output = filter(kernel, 1, data);

elseif method_idx == 2
    for i = 1:N
        output_f(:,i) = Fast_Smooth(data(:,i),windowWidth,3,1);
    end
    
    if T > 200
        temp_f = zeros(201,N);
        for i=1:N
            temp_f(:,i) = spline(linspace(0,1,T)',output_f(:,i),linspace(0,1,201)');
        end 
    end 
    
else
    warning('wrong index in Moving Average Routine.\n');
end

plot(grid,data(:,1),'-bo',grid,output_f(:,1),'-g*'); 
legend('data','output','location','best');