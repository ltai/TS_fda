%%% created on June, 30, 2016

function [ syn_para ] = Set_Synthetic_Data_Para

    syn_para.T = 201;    % number of time steps; Must be Odd; use 201!! 
    syn_para.N = 10;      % number of functions 
    syn_para.GA_choice = 3;  % group action choice; 2 is g(r)*dr; 
                             %                      3 is g(r)*sqrt(dr);
                          
    syn_para.noisy = 2;  % 0 is no noisy
                         % 1 is Gaussian noisy;   epsilon_i
                         % 2 is very noisy Gaussian; epsilon_i(t)
                       
    syn_para.normal_mu = 0;         % additive gaussian noise
    syn_para.normal_sigma = 0.2;                        

    %%%%%%%%%%%%%%%  set artificial data for g and h %%%%%%%%%%%%%%%%%    
    syn_para.generate_g = 65;   
                           % 10 Newton
                           % 15 Cos             '5*cos(10*pi*t)'
                           % 16                 'cos(10*pi*t)'
                           % 17                 '4*cos(10*pi*t)+5*cos(9*pi*t)'      
                           % 18                 '5*cos(9*pi*t)+10'
                           % 19                 'random coef w. Cos basis'
                           % 20 Sin             'sin(10*pi*t)'
                           % 21                 '5*sin(9*pi*t)'
                           % 22                 '5*sin(10*pi*t)+10' 
                           % 23                 '5*sin(9*pi*t)+10'
                           % 24                 'random coef w. Sin basis'
                           % 25 Modified F.S.   '5*cos(10*pi*t) + 5*sin(10*pi*t)'
                           % 26                 '4*cos(10*pi*t) + 6*sin(10*pi*t)'
                           % 27                 '5*cos(10*pi*t) + 5*sin(10*pi*t)+10'
                           % 28                 '4*cos(10*pi*t) + 6*sin(10*pi*t)+10'
                           % 29                 'random coef w. M.F.S. basis'
                           % 30 Legendre        'P_4','70*t.^4 - 140*t.^3 + 90*t.^2 - 20*t + 1'
                           % 31                 'P_7'
                           % 32                 'P_4 + P_7'
                           % 33                 '2*P_4 + 3*P_7'
                           % 34 
                           % 35                 random sLegendre poly, coef~exponential
                           % 40 Poly            random coef of poly. of degree n
                           % 41                 ' t.^8 '
                           % 42   
                           % 45                 'a poly but close to 5*cos(10*pi*t)'
                           % 46                 'a poly but close to cos(10*pi*t)'
                           % 50 Exponential     'exp(t)'
                           % 51       
                           % 55 Gaussian        '5*exp( -20*(t-0.5).^2 )'
                           % 56                 '3*exp( -20*(t-0.5).^2 )+3*exp( -20*(t+0.5).^2)'
                           % 65                 three peaks
                           % 66                 two peaks add up two exp fcns
                             
    syn_para.generate_h = 52;  
                           % 1 is no trend; h = 0
                           % 10 Newton
                           % 15 Cos             'cos(2*pi*t)'
                           % 16                 'cos(3*pi*t)'
                           % 17                 'cos(2*pi*t)+2*cos(pi*t)'      
                           % 18                 'cos(pi*t+pi/2)'
                           % 19                 'random coef w. Cos basis'
                           % 20 Sin             'sin(2*pi*t)'
                           % 21                 'sin(pi*t)'
                           % 22                 'sin(2*pi*t)+1' 
                           % 23                 'sin(pi*t)+1'
                           % 24                 'random coef w. Sin basis'
                           % 25 Modified F.S.   'cos(2*pi*t) + sin(2*pi*t)'
                           % 26                 '0.8*cos(2*pi*t) + 1.2*sin(2*pi*t)'
                           % 27                 'cos(2*pi*t) + sin(2*pi*t)+1'
                           % 28                 '0.8*cos(2*pi*t) + 1.2*sin(2*pi*t)+1'
                           % 29                 'random coef w. M.F.S. basis'
                           % 30 Legendre        '1'
                           % 31                 '2*t-1'
                           % 32                 '6*t.^2 -6*t +1'
                           % 33                 '20*t.^3 -30*t.^2 + 12*t - 1'
                           % 34                 '2*P_1 + 3*P_3'                 
                           % 35                 random sLegendre poly, coef~exponential
                           % 40 Poly            random coef of poly. of degree n
                           % 41                 '-0.1 + 0.3*t' linear trend
                           % 42                 '0.1-0.3*t+0.2*t.^2' quadratic
                           % 43                 ' 20(t-0.3)^2 *(t-0.7) ' cubic
                           % 44                 '(t-1).^2', slow increasing
                           % 45                 'a poly but very close to -log10(t+10^-3)'  
                           % 46                 '-2*t*(t-1)', parabola open down
                           % 47                 '(t-1).^2', slow decaying
                           % 48                 'a ploy but very close to exp(-5*t)'
                           % 49                 'a poly but very close to Gaussian'
                           % 50 Exponential     'exp(t)'
                           % 51                 '1.5*exp(-3*t)', decreasing
                           % 52                 '0.05*exp(3*t)-0.5', increasing
                           % 53                                
                           % 55 Gaussian        'exp( -20*(t-0.5).^2 )'
                           % 56                 '5*exp( -20*(t-0.5).^2 )'
                           % 60 logrithmic      '-log10(t + 10^-3)'
                           % 61                 '-2*log10(t + 10^-3)'

                            
    syn_para.generate_r = 2; % 1              smooth nice gamma
                             % 2              median gamma
                             % 3              ugly gamma
    %%%%%%%%%%%%%% No need to change below %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    syn_para.scaling = 0; % 0 is scaling 1; no update; always set to this. 
                          % 1 is scaling 1; but updating         %IC is 0.5
                          % 2 is exponantial noise; exprnd(1)    %IC is 0.5
                          % 3 is gaussian noise; N(1,0.1)        %IC is 0.5
                          % 4 is gaussian noice; N(2,0.1)        %IC is 0.5
                          % 5 is scaling 2; updating             %IC is 0.5
                          % 10 is scaling 2; no update;
                          
    syn_para.r_cond = 1;       % 0 is no cond.;
                               % 1 is KM{r^-1 _i} = r_id; (Default)
                               % 2 is AM{r^-1 _i} = r_id;
                      
    syn_para.g_cond = 0;     % 0 is no cond.; (Default)
                             % 1 is int = 0; 
                             % 2 is L2-norm = 1; 
                             % 3 is int = 0 and L2-norm = 1; 
                             % 4 is int = 1; 
                             % 5 is from invariance % forget the purpose
                             % 100 is int(g) = 0 and L2norm(g) = c1
                
    syn_para.h_cond = 0;    % 0 is no cond.;   (Default)
                            % 1 is int = 0; 
                            % 2 is L2-norm = 0.2; 
                            % 3 is int = 0 and L2-norm = 0.2; 
                            % 4 is L2-norm = 5; 
                            % 5 is int(abs) = 1;   
                            % 100 is int(h) = 0 and L2norm(h) = c2
                            % 101 is int(h) = c1 and L2norm(h) = c2; relate to f_i, g, and gam             
end

