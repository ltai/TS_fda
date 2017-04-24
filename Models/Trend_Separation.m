% this is the main function for trend separation, 
% continued on old version which is created on July, 27, 2014
% created on March, 21, 2015; all the regularization method be careful
%                             can be rewrited and the energy fcn associated
%                             with them need to be checked
% continued on July,1,2016; done :)

function [ est_data, final_est_h ] = Trend_Separation(para, true_data)

    tic;
    GA = para.GA_choice;
    max_iter = para.max_iter;
    T = true_data.T;
    N = true_data.N;
    f = true_data.f;
    t = true_data.t;
    
    est_h = zeros(T,max_iter);
    est_g = zeros(T,max_iter);
    gam = zeros(T,N,max_iter);
    gamI = zeros(T,N,max_iter);
    est_scalar = ones(max_iter,N);  % history purpose, scalar still in here
    
    Res = zeros(T,N);
    r_error = zeros(1,max_iter);
    L2cost = zeros(3,max_iter);
    g_error = zeros(1,max_iter);
    h_error = zeros(1,max_iter);
    ortho = zeros(1,max_iter);

%%%%%%%  generate initial estimation  %%%%%%
    %%%  true I.C. here %%%
    %   gam(:,:,1) = r;
    %   est_g(:,1) = g;
    %   est_h(:,1) = h;
    
    %%%%%% initialized gam %%%%%
    for i = 1:N
        gam(:,i,1) = t;
    end
    
    %%%%%% initialized g %%%%%
    % est_g(:,1) = mean(true_data.f,2); % turns out to be bad
    % est_g(:,1) = f(:,5);              % if not putting KM of gam
    temp_L2norm = zeros(1,N);           % to the initialization of g...
    mean_f = mean(true_data.f,2);       % Feb. 3, 2016
    for i = 1:N
        temp_L2norm(i) = L2norm(t, f(:,i) - mean_f);
    end
    [~,min_idx] = min(temp_L2norm);  fprintf(' initial g idx: %d',min_idx);     
    est_g(:,1) = f(:,min_idx);
    
    %%%%%% initialized h %%%%%
    % already set to be zero
    
%%%%%%%%%%%  run Coordinate Descent  %%%%%%%%%%%%%%
%%%%%%%%%%%  run Coordinate Descent  %%%%%%%%%%%%%%
    
    for j = 2:max_iter
        % keyboard;
        fprintf('\n %d-th iter: ',j);
        
%%%%%%%%%%    update r    %%%%%%%%%%%%
%%%%%%%%%%    update r    %%%%%%%%%%%%
        if para.update_r == 0
            if true_data.use_data == 0
                disp('wrong parameter setting'); break;        
            elseif true_data.use_data == 1
                gam(:,:,j) = true_data.r;
            end
                
        elseif para.update_r == 2 || para.update_r == 3
            lambda = 0;
            for i = 1:N
                if para.update_r == 2
                    gam(:,i,j) = MyDynamicProgramming( est_g(:,j-1)', (f(:,i) - est_h(:,j-1))'/est_scalar(j-1,i),lambda,2)';
                elseif para.update_r == 3
                    gam(:,i,j) = MyDynamicProgramming( est_g(:,j-1)', (f(:,i) - est_h(:,j-1))'/est_scalar(j-1,i),lambda,3)';
                end
            end
         
        % KM need to impose at each iteration, only last, very first iteration will
        % fail, by Feb. 2, 2016 testing result
        % converging case: KM at each iteration makes the cost function
        % down at each iteration. 
            if para.est_r_cond == 0
            elseif para.est_r_cond == 1
                gam(:,:,j) = Gamma_with_inverse_kmean_id(gam(:,:,j));
            elseif para.est_r_cond == 2
                gam(:,:,j) = Gamma_with_inverse_Amean_id(gam(:,:,j));
            else
                warning('wrong index of r_cond');
            end   
        else
            warning('wrong index in para.update_r');
        end
        fprintf('done r; ');       
        
%%%%%%%%%%%%    update g, imediate update   %%%%%%%%%%%%
%%%%%%%%%%%%    update g, imediate update   %%%%%%%%%%%%
        if para.update_g == 0
            est_g(:,j) = true_data.g;
            
        elseif para.update_g == 1 % non-para average
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                est_g(:,j) = est_g(:,j) + ...
                        MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
            end
            est_g(:,j) = est_g(:,j)/N;
            
        elseif  para.update_g == 15 % para average, Cos;
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                Res(:,i) = MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
                est_g(:,j) = est_g(:,j) + L2_Cos_Approximation(Res(:,i),para.g_start_idx,para.g_end_idx); 
            end
            est_g(:,j) = est_g(:,j)/N;
                
        elseif  para.update_g == 20 % para average, Sin
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                Res(:,i) = MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
                est_g(:,j) = est_g(:,j) + L2_Sin_Approximation(Res(:,i),para.g_start_idx,para.g_end_idx);
            end
            est_g(:,j) = est_g(:,j)/N;
            
        elseif para.update_g == 25 % para average, Modified Fourier Basis
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                Res(:,i) = MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
                est_g(:,j) = est_g(:,j) + L2_MFS_Approximation(Res(:,i),para.g_start_idx,para.g_end_idx);
            end
            est_g(:,j) = est_g(:,j)/N;
            
        elseif  para.update_g == 30 % para average, Shifted Legendre Basis   
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                Res(:,i) = MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
                est_g(:,j) = est_g(:,j) + ...                          
                    L2_Shift_Legendre_Approximation(Res(:,i),para.g_start_idx,para.g_end_idx);
                    %L2_Shift_Legendre_Approximation(Res(:,i),5,20);
            end
            est_g(:,j) = est_g(:,j)/N;
            
        elseif para.update_g == 40 % nonpara average with ortho. to para. h
            for i = 1:N
                gamI(:,i,j) = invertGamma(gam(:,i,j));
                est_g(:,j) = est_g(:,j) + ...
                        MyGroupAction(t, (f(:,i) - est_h(:,j-1))/est_scalar(j-1,i), gamI(:,i,j),GA);
            end
            est_g(:,j) = est_g(:,j)/N;
            
            proj_g = zeros(T,1);
            if para.update_h == 15
                proj_g = L2_Cos_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
           
            elseif para.update_h == 20
                proj_g = L2_Sin_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
                
            elseif para.update_h == 25
                proj_g = L2_MFS_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
                
            elseif para.update_h == 30
                proj_g = L2_Shift_Legendre_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
                
            else
                warning('wrong in nonpara avg g with ortho. to para h\n');
            end
            est_g(:,j) = est_g(:,j) - proj_g;      
        end
        fprintf('done g; ');    
        
%%%%%%%%%%%    update h,  imetidiate update !     %%%%%%%%%%%%
%%%%%%%%%%%    update h,  imetidiate update !     %%%%%%%%%%%%
        if para.update_h == 0   % no update h
            est_h(:,j) = true_data.h;
            
        elseif para.update_h == 1 %%% avegage method; no regularization 
            for i = 1:N
                est_h(:,j) = est_h(:,j) + f(:,i) ...
                    - est_scalar(j-1,i)*MyGroupAction(t, est_g(:,j), gam(:,i,j),GA);
            end
            est_h(:,j) = est_h(:,j)/N;
                   
        elseif para.update_h == 15   %%%% para average, Cos Basis;    
            for i = 1:N
                warp_g = est_scalar(j-1,i)*MyGroupAction(t,est_g(:,j),gam(:,i,j),GA);
                Res(:,i) = f(:,i) - warp_g;                            
                est_h(:,j) = est_h(:,j) + L2_Cos_Approximation(Res(:,i),para.h_start_idx,para.h_end_idx);
            end
            est_h(:,j) = est_h(:,j)/N;
                 
        elseif para.update_h == 20   %%%% para average, Sin Basis
            for i = 1:N
                warp_g = est_scalar(j-1,i)*MyGroupAction(t,est_g(:,j),gam(:,i,j),GA);
                Res(:,i) = f(:,i) - warp_g;                               
                est_h(:,j) = est_h(:,j) + L2_Sin_Approximation(Res(:,i),para.h_start_idx,para.h_end_idx);
            end
            est_h(:,j) = est_h(:,j)/N;
            
        elseif para.update_h == 25   %%%% para average, Modified Fourier Basis    
             for i = 1:N
                warp_g = est_scalar(j-1,i)*MyGroupAction(t,est_g(:,j),gam(:,i,j),GA);
                Res(:,i) = f(:,i) - warp_g;                               
                est_h(:,j) = est_h(:,j) + L2_MFS_Approximation(Res(:,i),para.h_start_idx,para.h_end_idx);
            end
            est_h(:,j) = est_h(:,j)/N;
            
        elseif para.update_h == 30 % para average, Shifted Legendre Basis   
            for i = 1:N
                warp_g = est_scalar(j-1,i)*MyGroupAction(t,est_g(:,j),gam(:,i,j),GA);
                Res(:,i) = f(:,i) - warp_g;                          
                est_h(:,j) = est_h(:,j) + ...  
                    L2_Shift_Legendre_Approximation(Res(:,i),para.h_start_idx,para.h_end_idx);
            end
            est_h(:,j) = est_h(:,j)/N;
            
        elseif para.update_h == 40 % non-para average with ortho. to para.g
            for i = 1:N
                est_h(:,j) = est_h(:,j) + f(:,i) ...
                    - est_scalar(j-1,i)*MyGroupAction(t, est_g(:,j), gam(:,i,j),GA);
            end
            est_h(:,j) = est_h(:,j)/N;
            
            proj_h = zeros(T,1);
            if para.update_g == 15
                proj_h = L2_Cos_Approximation(est_h(:,j),para.g_start_idx,para.g_end_idx);
                
            elseif para.update_g == 20
                proj_h = L2_Sin_Approximation(est_h(:,j),para.g_start_idx,para.g_end_idx);
                
            elseif para.update_g == 25
                proj_h = L2_MFS_Approximation(est_h(:,j),para.g_start_idx,para.g_end_idx);
                
            elseif para.update_g == 30
                proj_h = L2_Shift_Legendre_Approximation(est_h(:,j),para.g_start_idx,para.g_end_idx);
                
            else
                warning('wrong in nonpara avg g with ortho. to para h\n');
            end
            est_h(:,j) = est_h(:,j) - proj_h;
        end
        fprintf('done h; ');           
    end % end of Coordinate Descent
    timecost = toc;

%%%%%%%%%% compute error %%%%%%%%%%%%%%%%% 
%%%%%%%%%% compute error %%%%%%%%%%%%%%%%%
    fprintf('\n Computing Error');
    for j = 1:max_iter
        fprintf('\n %d-th iter: ',j);
        if true_data.use_data == 1 % synthetic data has ground truth
            g_error(j) = L2norm(t, est_g(:,j) - true_data.g)/L2norm(t,true_data.g);
            h_error(j) = L2norm(t, est_h(:,j) - true_data.h)/L2norm(t,true_data.h);
            r_error(j) = 0;
            for i = 1:N
                r_error(j) = r_error(j) + L2norm(t,gam(:,i,j) - true_data.r(:,i));
            end
            r_error(j) = r_error(j)/N;
            true_ortho = Check_Orthogonality(true_data.g,true_data.h,true_data.t);
            ortho(j) = Check_Orthogonality(est_g(:,j),est_h(:,j), true_data.t);
        else
            ortho(j) = Check_Orthogonality(est_g(:,j),est_h(:,j), linspace(0,1,true_data.T));
        end
        
        if j == 1
            L2cost(:,1) = Compute_Cost(para, f, est_g(:,1), gam(:,:,1), est_h(:,1),est_scalar(1,:));
        else
            % for L2 error of updating gam
            L2cost(1,j) = Compute_Cost(para, f, est_g(:,j-1), gam(:,:,j), est_h(:,j-1),est_scalar(j-1,:));
            if L2cost(1,j) - L2cost(3,j-1) > 10/T
                fprintf('increase; ');
            else
                fprintf('decrease; ');
            end
            
            % for L2 error of updating g
            L2cost(2,j) = Compute_Cost(para, f, est_g(:,j), gam(:,:,j), est_h(:,j-1),est_scalar(j-1,:));
            if L2cost(2,j) - L2cost(1,j) > 10/T
                fprintf('increase; ');
            else
                fprintf('decrease; ');
            end
            
            % for L2 error of updating h
            L2cost(3,j) = Compute_Cost(para, f, est_g(:,j), gam(:,:,j), est_h(:,j),est_scalar(j-1,:));
            if L2cost(3,j) - L2cost(2,j) > 10/T
                fprintf('increase; ');
            else
                fprintf('decrease; ');
            end
        end
    end    
    
    if true_data.use_data == 1
        scalar = ones(1,N);
        syn_L2cost = Compute_Cost(para, f, true_data.g, true_data.r, true_data.h, scalar);
    end

    %%%%%  output estimated data  %%%%%%
    if true_data.use_data == 1
        est_data.syn_L2cost = syn_L2cost;
        est_data.g_error = g_error;
        est_data.h_error = h_error;
        est_data.r_error = r_error;
        est_data.true_ortho = true_ortho;
    end
    
    est_data.gam = gam;
    est_data.gamI = gamI;
    est_data.est_g = est_g;
    est_data.est_h = est_h;
    
    est_data.timecost = timecost;
    est_data.ortho = ortho;
    est_data.L2cost = L2cost;
    est_data.final_L2cost = L2cost(end,end);
    final_est_h = est_h(:,end);
end
