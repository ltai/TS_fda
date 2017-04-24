% created on Jan. 7, 2016
% finished on Jan. 8, 2016! successful tested.

function [est_data] = Simple_Separation(para, true_data)

    %%%% initialization
    T = true_data.T;
    N = true_data.N;
    f = true_data.f;
    max_iter = para.max_iter;
    t = true_data.t;
    
    h_start_idx = para.h_start_idx;
    h_end_idx = para.h_end_idx;

    est_g = zeros(T,max_iter);
    est_h = zeros(T,max_iter);
    L2cost = zeros(max_iter,1);
  
    est_g(:,1) = mean(f,2);

    for j = 2:max_iter      % this model don't need to do iterations
        % update est g
        temp_vec = zeros(T,1);
        for i = 1:N
            temp_vec = temp_vec + f(:,i) - est_h(:,j-1);
        end
        est_g(:,j) = temp_vec/N;
        
        if para.update_h == 15
            est_g(:,j) = est_g(:,j) - L2_Cos_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
        elseif para.update_h == 20
            est_g(:,j) = est_g(:,j) - L2_Sin_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
        elseif para.update_h == 25
            est_g(:,j) = est_g(:,j) - L2_MFS_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
        elseif para.update_h == 30
            est_g(:,j) = est_g(:,j) - L2_Shift_Legendre_Approximation(est_g(:,j),para.h_start_idx,para.h_end_idx);
        end
        
        % update est h
        temp_vec = zeros(T,1);
        for i = 1:N
            if para.update_h == 15
                temp_vec = temp_vec + L2_Cos_Approximation(f(:,i), h_start_idx, h_end_idx);
            elseif para.update_h == 20 
                temp_vec = temp_vec + L2_Sin_Approximation(f(:,i), h_start_idx, h_end_idx);
            elseif para.update_h == 25
                temp_vec = temp_vec + L2_MFS_Approximation(f(:,i), h_start_idx, h_end_idx);
            elseif para.update_h == 30
                temp_vec = temp_vec + L2_Shift_Legendre_Approximation(f(:,i), h_start_idx, h_end_idx);
            end
        end
        est_h(:,j) = temp_vec/N;
    end

    %%% compute error
    for j = 1:max_iter
        temp_L2cost = 0;
        for i = 1:N
            temp_L2cost = temp_L2cost + SquareL2norm(t,f(:,i) - est_g(:,j) - est_h(:,j) );
        end
        L2cost(j) = temp_L2cost/N;
    end
    
    %%% compute error if synthetic data
    if true_data.use_data == 1
        g_error = L2norm(t, est_g(:,end) - true_data.g)/L2norm(t,true_data.g);   
        h_error = L2norm(t, est_h(:,end) - true_data.h)/L2norm(t,true_data.h);
        ortho = Check_Orthogonality(est_g(:,end),est_h(:,end), true_data.t);
    else
        ortho = Check_Orthogonality(est_g(:,end),est_h(:,end), linspace(0,1,true_data.T));
    end

    est_data.est_g = est_g;
    est_data.est_h = est_h;
    est_data.L2cost = L2cost;
    est_data.final_L2cost = L2cost(end);
    if true_data.use_data == 1
        est_data.g_error = g_error;
        est_data.h_error = h_error;
    end
    est_data.ortho = ortho;
end

