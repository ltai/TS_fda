%%%% created on June, 30,2016, successful

function [ true_data ] = Import_Data(data)

    smooth_choice = 0;    % 0 is no smooth (Default); 1 is smooth

if nargin == 1
    f = data;
    % no smoothing
    [T,N] = size(f);
    t = linspace(0,1,T)';
    
elseif nargin == 0
%%%%%%%%%%%%%  set real or artificial data   %%%%%%%%%%%%%%%%%%%%%%%%%%                     
   real_data = 61;       % 1 is male growth velocity first 10 
                        % 2 is male growth velocity; whole data    
                        % 3 is 
                        % 4 is
                        % 5 is river_flow 1915 - 1924
                        % 6 is 
                        % 7 is 
                        % 8 is 
                        % 9 is 
                        % 10 is Spectrometry Data
                        % 20 is 
                        % 30 is U.S. bank Stock price smooth
                        % 31 is US_Bank_Stock_Price_Difference_smooth
                        % 35 is Florida Monthly Temparature
                        % 37 is Florida daily average temparature
                        % 40 is Florida Monthly Sea Level (mm),
                        % 41 is Florida Monthly Sea level, (mm), smooth 13
                        % 42 is Florida Monthly Sea Level (cm),
                        % 43 is Florida Monthly Sea level, (cm), smooth 12
                        % 50 is Tallahassee Daily Temperature 
                        % 51 is Tallahassee Daily Temperature, smooth 31
                        % 52 is Tallahassee Daily Temperature, smooth 11
                        % 55 is Tallahassee Daily Precipitation
                        % 56 is Tallahassee Daily Precipitation box smooth
                        % 57 is Tallahassee Daily Precipitation gauss smooth
                        % 60 is currency exchange rate, difference rate, to USD
                        % 61 is Currency_Exchange_To_USD_2015_10to12_smooth
                        % 63 is USD_Euro_exchange_spline_smooth
                        % 65 is U.S. electrical price, 2009 - 2014,
                        % 66 is U.S. electrical price difference, 2005 - 2010, box smooth
                        % 67 is U.S. electrical price difference, 2005 - 2010, gauss smooth
                        % 69 is US_Electricity_Price_2009_2014_remove_outlier_6
                        % 100 artificial g=5cos(10pi*t),h=-log10(t+0.001),eps~N(0,1) 

        if real_data == 1
            load growth_male_vel_10;
            
        elseif real_data == 2
            load growth_male_vel;
            
        elseif real_data == 5
            load river_flow_1915_1924; 
            
        elseif real_data == 10
            load Spectrometry_Data_First20;
            
        elseif real_data == 20
            load accelsig_U1;
            
        elseif real_data == 30
            load Bank_Stock_Price_smooth;
        
        elseif real_data == 31
            load US_Bank_Stock_Price_Difference_smooth;
            
        elseif real_data == 35
            load Florida_Temparature;
            
        elseif real_data == 37
            load TLH_AverageTemperature_smooth_10;
            
        elseif real_data == 40
            load Florida_Sea_Level_mm_2005_2014;
            
        elseif real_data == 41
            load Florida_Sea_Level_mm_2005_2014_smooth;
            
        elseif real_data == 42
            load Florida_Sea_Level_cm_2005_2014;
            
        elseif real_data == 43
            load Florida_Sea_Level_cm_2005_2014_smooth;
            
        elseif real_data == 50
            load Tallahassee_Daily_Temperature;
            
        elseif real_data == 51
            load Tallahassee_Daily_Temperature_smooth31;
            
        elseif real_data == 52
            load Tallahassee_Daily_Temperature_smooth11;    
            
        elseif real_data == 55
            load Tallahassee_Daily_Precipitation_smooth;
            
        elseif real_data == 56
            load Tallahassee_Daily_Precipitation_box_smooth;
            
        elseif real_data == 57
            load Tallahassee_Daily_Precipitation_gauss_smooth;
            
        elseif real_data == 60    
            load Currency_Exchange_Difference_To_USD_2014_smooth;
            
        elseif real_data == 61
            load Currency_Exchange_To_USD_2015_10to12_smooth;
            
        elseif real_data == 63 
            load USD_Euro_exchange_spline_smooth;
            
        elseif real_data == 65
            load US_Electricity_Price_2009_2014;
            
        elseif real_data == 66
            load US_Electricity_Price_Difference_2005_2010_box_smooth;
            
        elseif real_data == 67
            load US_Electricity_Price_Difference_2005_2010_gauss_smooth;
            
        elseif real_data == 68
            load US_Electricity_Price_2009_2014_remove_outlier_7;
            
        elseif real_data == 69
            %load US_Electricity_Price_2009_2014_remove_outlier_5;
            load US_Electricity_Price_2009_2014_remove_outlier_6;
            
        elseif real_data == 100
            load Artificial_fi_cos_log;
            
        else
            warning('wrong index when selecting real data file\n');
        end
        
        % smoothing 
        if smooth_choice == 0
            [T,N] = size(f);
            t = linspace(0,1,T)';
            
        elseif smooth_choice == 1
            [temp_T,N] = size(f);
            temp_t = linspace(0,1,temp_T)';
            T = 121; % need odd points.
            t = linspace(0,1,T)';
            temp_f = zeros(T,N);
            for i = 1:N
                temp_f(:,i) = spline(temp_t,f(:,i),t);
            end
            clear f;
            f = temp_f;
        end
end


       
        %%%%%%%%% output data
        true_data.T = T;
        true_data.N = N;
        true_data.use_data = 0;     % 0 is real or load data 
    
        true_data.t = t;
        true_data.f = f;
        mean_f = mean(f,2);
        true_data.mean_f = mean_f;
end

