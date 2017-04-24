% this is the main function for trend separation, 
% continued on old version which is created on July, 27, 2014
% created on March, 21, 2015

function [ true_data ] = Generate_Synthetic_Data( syn_para )
        
    GA = syn_para.GA_choice;

    T = syn_para.T;
    N = syn_para.N;
    t = linspace(0,1,T)';      % time step, t

%%%%%%%%%%%%%%%%%%% generate gamma %%%%%%%%%%%%%%%%%%%%%%%%%
    
    if syn_para.generate_r == 1
        r = get_gam(T,N,4);
    elseif syn_para.generate_r == 2
        r = get_gam(T,N,5);
    elseif syn_para.generate_r == 3
        r = get_gam(T,N,6);
    else
        warning('wrong syn_para.generate_h index'); 
    end
    
    if syn_para.r_cond == 0
    elseif syn_para.r_cond == 1
        r = Gamma_with_inverse_kmean_id(r);
    elseif syn_para.r_cond == 2
        r = Gamma_with_inverse_Amean_id(r);
    else
        warning('wrong idx selection para.r_cond');
    end

%%%%%%%%%%%%%%%%  generate g %%%%%%%%%%%%%%%%%%%%%%%
    g(:,1) = zeros(T,1);

    if syn_para.generate_g == 0  % won't be this case,
    elseif syn_para.generate_g == 10 % Newton, but no longer needed
    elseif syn_para.generate_g == 15
        fcn = inline('5*cos(10*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 16
        fcn = inline('cos(10*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 17
        fcn = inline('4*cos(10*pi*t)+5*cos(9*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 18
        fcn = inline('5*cos(9*pi*t)+10','t'); g = fcn(t);

    elseif syn_para.generate_g == 19
        rand_coef = 0.1 + randn(100,2);
        g = zeros(T,1);

        for i = syn_para.g_start_idx : syn_para.g_end_idx
            g = g + rand_coef(i,1)*cos(pi*(i-1)*t);
        end

    elseif syn_para.generate_g == 20
        fcn = inline('sin(10*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 21
        fcn = inline('5*sin(9*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 22
        fcn = inline('5*sin(10*pi*t)+10','t'); g = fcn(t);

    elseif syn_para.generate_g == 23
        fcn = inline('5*sin(9*pi*t)+10','t'); g = fcn(t);

    elseif syn_para.generate_g == 24
        rand_coef = 0.1 + randn(100,2);
        g = zeros(T,1);

        for i = syn_para.g_start_idx : syn_para.g_end_idx
            g = g + rand_coef(i,1)*sin(pi*(i-1)*t);
        end

    elseif syn_para.generate_g == 25
        fcn = inline('5*cos(10*pi*t) + 5*sin(10*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 26
        fcn = inline('4*cos(10*pi*t) + 6*sin(10*pi*t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 27
        fcn = inline('5*cos(10*pi*t) + 5*sin(10*pi*t)+10','t'); g = fcn(t);

    elseif syn_para.generate_g == 28
        fcn = inline('4*cos(10*pi*t) + 6*sin(10*pi*t)+10','t'); g = fcn(t);

    elseif syn_para.generate_g == 29
        rand_coef = 0.1 + randn(100,2);
        g = zeros(T,1);

        for i = syn_para.g_start_idx : syn_para.g_end_idx
            g = g + rand_coef(i,1)*cos(2*pi*(i-1)*t) ...
                + rand_coef(i,2)*sin(2*pi*(i-1)*t);
        end

    elseif syn_para.generate_g == 30
        fcn = inline('70*t.^4 - 140*t.^3 + 90*t.^2 - 20*t + 1','t'); g = fcn(t);

    elseif syn_para.generate_g == 31
        fcn = inline('3432*t.^7 - 12012*t.^6 + 16632*t.^5 - 11550*t.^4 + 4200*t.^3 - 756*t.^2 + 56*t - 1','t');
        g = fcn(t);

    elseif syn_para.generate_g == 32
        fcn1 = inline('70*t.^4 - 140*t.^3 + 90*t.^2 - 20*t + 1','t');
        fcn2 = inline('3432*t.^7 - 12012*t.^6 + 16632*t.^5 - 11550*t.^4 + 4200*t.^3 - 756*t.^2 + 56*t - 1','t');
        g = fcn1(t) + fcn2(t);

    elseif syn_para.generate_g == 33
        fcn1 = inline('70*t.^4 - 140*t.^3 + 90*t.^2 - 20*t + 1','t');
        fcn2 = inline('3432*t.^7 - 12012*t.^6 + 16632*t.^5 - 11550*t.^4 + 4200*t.^3 - 756*t.^2 + 56*t - 1','t');
        g = 2*fcn1(t) + 3*fcn2(t);

    elseif syn_para.generate_g == 35
        temp = load('Legendre_Poly_1001','Legendre_Poly');
        Legendre_Poly = temp.Legendre_Poly;
        [temp_T,~] = size(Legendre_Poly);
        tt = linspace(0,1,temp_T)';

        temp_g = zeros(temp_T,1);
        rand_coef = 0.05 + randn(20,1);

        for i = syn_para.g_start_idx : syn_para.g_end_idx
            temp_g = temp_g + rand_coef(i)*Legendre_Poly(:,i);
        end

        t = linspace(0,1,T)';
        g = spline(tt,temp_g,t);

    elseif syn_para.generate_g == 40
        Degree = syn_para.Degree;
        a = -1/2 + 1*rand(Degree+1,1);
        for k = 1:Degree+1
            g(:,1) = g(:,1) + a(k,1)*t.^(k-1);
        end

    elseif syn_para.generate_g == 41
        fcn = inline('t.^8 ','t'); g = fcn(t);

    elseif syn_para.generate_g == 45
        load g_ploy_close_to_5cos10pit;

    elseif syn_para.generate_g == 46
        load g_ploy_close_to_cos10pit;

    elseif syn_para.generate_g == 50
        fcn = inline('exp(t)','t'); g = fcn(t);

    elseif syn_para.generate_g == 55
        fcn = inline('5*exp( -20*(t-0.5).^2 )','t'); g = fcn(t);

    elseif syn_para.generate_g == 56
        fcn = inline('3*exp( -15*(t-0.7).^2 )+ 1*exp( -15*(t+0.8).^2)','t');
        g = fcn(t);
    elseif syn_para.generate_g == 65
        g = 5*sin(5*pi*t).*(0.25-(t-0.5).^2); 
        
    elseif syn_para.generate_g == 66
        g = 2*exp(-0.8*(10*t-7.5).^2) + 2*exp(-0.8*(10*t-2.5).^2)-0.85;
    else
        warning('wrong index in generating g\n');
    end

%%%%%%%%%%%%%%%%%%%%% generate h %%%%%%%%%%%%%%%%%%%%%%%%%
    h(:,1) = zeros(T,1);

    if syn_para.generate_h == 0
    elseif syn_para.generate_h == 1
        h = 0*t;

    elseif syn_para.generate_h == 10

    elseif syn_para.generate_h == 15
        trend = inline('cos(2*pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 16
        trend = inline('cos(3*pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 17
        trend = inline('cos(2*pi*t)+2*cos(pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 18
        trend = inline('cos(pi*t+pi/2)','t'); h = trend(t);

    elseif syn_para.generate_h == 19
        rand_coef = 0.1 + randn(100,2);
        h = zeros(T,1);

        for i = syn_para.h_start_idx : syn_para.h_end_idx
            h = h + rand_coef(i,1)*cos(pi*(i-1)*t);
        end

    elseif syn_para.generate_h == 20
        trend = inline('sin(2*pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 21
        trend = inline('sin(pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 22
        trend = inline('sin(2*pi*t)+1','t'); h = trend(t);

    elseif syn_para.generate_h == 23
        trend = inline('sin(pi*t)+1','t'); h = trend(t);

    elseif syn_para.generate_h == 24
        rand_coef = 0.1 + randn(100,2);
        h = zeros(T,1);
        for i = syn_para.h_start_idx : syn_para.h_end_idx
            h = h + rand_coef(i,1)*sin(pi*(i-1)*t);
        end

    elseif syn_para.generate_h == 25
        trend = inline('cos(2*pi*t) + sin(2*pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 26
        trend = inline('0.8*cos(2*pi*t) + 1.2*sin(2*pi*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 27
        trend = inline('cos(2*pi*t) + sin(2*pi*t)+1','t'); h = trend(t);

    elseif syn_para.generate_h == 28
        trend = inline('0.8*cos(2*pi*t) + 1.2*sin(2*pi*t)+1','t'); h = trend(t);

    elseif syn_para.generate_h == 29
        rand_coef = 0.1 + randn(100,2);
        h = zeros(T,1);
        for i = syn_para.h_start_idx : syn_para.h_end_idx
            h = h + rand_coef(i,1)*cos(2*pi*(i-1)*t) ...
                + rand_coef(i,2)*sin(2*pi*(i-1)*t);
        end

    elseif syn_para.generate_h == 30
        trend = inline('0*t+1','t'); h = trend(t);

    elseif syn_para.generate_h == 31
        trend = inline('2*t-1','t'); h = trend(t);

    elseif syn_para.generate_h == 32
        trend = inline('6*t.^2 -6*t +1','t'); h = trend(t);

    elseif syn_para.generate_h == 33
        trend = inline('20*t.^3 -30*t.^2 + 12*t - 1','t'); h = trend(t);

    elseif syn_para.generate_h == 34
        fcn1 = inline('2*t-1','t');
        fcn2 = inline('20*t.^3 -30*t.^2 + 12*t - 1','t');
        h = 2*fcn1(t) + 3*fcn2(t);

    elseif syn_para.generate_h == 35
        temp = load('Legendre_Poly_1001','Legendre_Poly');
        Legendre_Poly = temp.Legendre_Poly;
        [temp_T,~] = size(Legendre_Poly);
        tt = linspace(0,1,temp_T)';

        temp_h = zeros(temp_T,1);
        rand_coef = 0.05 + randn(20,1);

        for i = syn_para.h_start_idx : syn_para.h_end_idx
            temp_h = temp_h + rand_coef(i)*Legendre_Poly(:,i);
        end

        t = linspace(0,1,T)';
        h = spline(tt,temp_h,t);

    elseif syn_para.generate_h == 40
        Degree = syn_para.Degree;
        a = -1/2 + 1*rand(Degree+1,1);
        for k = 1:Degree+1
            h(:,1) = h(:,1) + a(k,1)*t.^(k-1);
        end

    elseif syn_para.generate_h == 41
        trend = inline('-0.1 + 0.3*t','t'); h = trend(t);

    elseif syn_para.generate_h == 42
        trend = inline('0.1 - 0.3*t + 0.2*t.^2','t'); h = trend(t);

    elseif syn_para.generate_h == 43
        trend = inline('20*(t-0.3).^2 .*(t-0.7)','t'); h = trend(t);

    elseif syn_para.generate_h == 44
        trend = inline('-(t-1).^2','t'); h = trend(t);

    elseif syn_para.generate_h == 45
        load h_ploy_close_to_negative_log;

    elseif syn_para.generate_h == 46
        trend = inline('-2*t.*(t-1)','t'); h = trend(t);

    elseif syn_para.generate_h == 47
        trend = inline('(t-1).^2','t'); h = trend(t);

    elseif syn_para.generate_h == 48
        load h_poly_close_to_exp_negative5t;

    elseif syn_para.generate_h == 49
        load h_poly_close_to_gaussian;

    elseif syn_para.generate_h == 50
        trend = inline('exp(t)','t'); h = trend(t);

    elseif syn_para.generate_h == 51
        trend = inline('1.5*exp(-3*t)','t'); h = trend(t);

    elseif syn_para.generate_h == 52
        trend = inline('0.05*exp(3*t)-0.5','t'); h = trend(t);

    elseif syn_para.generate_h == 55
        trend = inline('exp( -20*(t-0.5).^2 )','t'); h = trend(t);

    elseif syn_para.generate_h == 56
        trend = inline('5*exp( -20*(t-0.5).^2 )','t'); h = trend(t);

    elseif syn_para.generate_h == 60
        trend = inline('-log10(t + 10^-3)','t'); h = trend(t);

    elseif syn_para.generate_h == 61
        trend = inline('-2*log10(t + 10^-3)','t'); h = trend(t);

    else
        fprintf('wrong index in generating h\n');
    end

    %%%% add scalar noice to g, no longer needed %%%%
    scalar = ones(1,N);

    if syn_para.scaling == 2  % exponential scalar
        scalar = exprnd(1, [1 N]);

    elseif syn_para.scaling == 3 % Gaussian scalar
        scalar = 1 + 0.1*randn(1,N);

    elseif syn_para.scaling == 4
        scalar = 2 + 0.1*randn(1,N);

    elseif syn_para.scaling == 5
        scalar = ones(1,N) + 1;

    elseif syn_para.scaling == 10
        scalar = zeros(1,N) + 2;
    end

%%%%%%%%%  additive Gaussian noise %%%
    rng(1);
    error = syn_para.normal_mu + syn_para.normal_sigma*randn(1,N);
    rng(1);
    error2 = syn_para.normal_mu + syn_para.normal_sigma*randn(T,N);

%%%%%%%%%%%%  generate observations f_i(t)  %%%%%%%%%
    f = zeros(T,N);
    for i = 1:N
        if syn_para.noisy == 1
            f(:,i) = scalar(i)*MyGroupAction(t, g, r(:,i),GA) + h(:,1) + error(1,i);
        elseif syn_para.noisy == 2
            f(:,i) = scalar(i)*MyGroupAction(t, g, r(:,i),GA) + h(:,1) + error2(:,i);
        else
            f(:,i) = scalar(i)*MyGroupAction(t, g, r(:,i),GA) + h(:,1);
        end
    end

%%%%  passing generated parameters to true data  %%%
    true_data.T = T;
    true_data.N = N;
    true_data.scalar = scalar;  
    true_data.use_data = 1;     % 1 is synthetic data 
    
    true_data.t = t;
    true_data.g = g;
    true_data.h = h;
    true_data.r = r;
    
    true_data.f = f;
    mean_f = mean(f,2);
    true_data.mean_f = mean_f; 
end