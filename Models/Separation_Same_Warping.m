% created on Jan. 7, 2016
% continued on Jan. 12, 2016, during the process

function [est_data] = Separation_Same_Warping(para, true_data)

    T = true_data.T;
    N = true_data.N;
    f = true_data.f;
    max_iter = para.max_iter;
    t = true_data.t;
    
    h_start_idx = para.h_start_idx;
    h_end_idx = para.h_end_idx;
    
    est_q_g_plus_h = zeros(T,max_iter);
    est_g_plus_h = zeros(T,max_iter);
    est_g = zeros(T,max_iter);
    est_h = zeros(T,max_iter);
    gam = zeros(T,N,max_iter);
    gamI = zeros(T,N,max_iter);
    q_f = zeros(T,N);
    L2cost = zeros(max_iter,1);
    SRVF_L2cost = zeros(max_iter,1);
    
    for i = 1:N
        q_f(:,i) = SRVF(t,f(:,i));
    end
    est_q_g_plus_h(:,1) = mean(q_f,2);
    
    for i = 1:N
        gam(:,i,1)  = t;
    end
    
    for j = 2:max_iter        
        
        %%%%% update gam
        lambda = 0;
        for i = 1:N
            gam(:,i,j) = MyDynamicProgramming( est_q_g_plus_h(:,j-1)', q_f(:,i)',lambda,3)';
        end
        gam(:,:,j) = Gamma_with_inverse_kmean_id(gam(:,:,j));
        
        %%%%% update q_g_plus_h        
        for i = 1:N
            gamI(:,i,j) = invertGamma(gam(:,i,j));     
            est_q_g_plus_h(:,j) = est_q_g_plus_h(:,j) + MyGroupAction(t, q_f(:,i), gamI(:,i,j),3);                
        end
        est_q_g_plus_h(:,j) = est_q_g_plus_h(:,j)/N;    
        
        est_g_plus_h(:,j) = InverseSRVF(t,est_q_g_plus_h(:,j), mean(f(1,:)));     
        
        %%%%% do projection into est_g and est_h
        for i = 1:N
            if para.update_h == 15
                est_h(:,j) = est_h(:,j) + L2_Cos_Approximation(est_g_plus_h(:,j), h_start_idx, h_end_idx);
            elseif para.update_h == 20
                est_h(:,j) = est_h(:,j) + L2_Sin_Approximation(est_g_plus_h(:,j), h_start_idx, h_end_idx);
            elseif para.update_h == 25
                est_h(:,j) = est_h(:,j) + L2_MFS_Approximation(est_g_plus_h(:,j), h_start_idx, h_end_idx);
            elseif para.update_h == 30
                est_h(:,j) = est_h(:,j) + L2_Shift_Legendre_Approximation(est_g_plus_h(:,j), h_start_idx, h_end_idx);
            end
        end
        est_h(:,j) = est_h(:,j)/N;
        
        est_g(:,j) = est_g_plus_h(:,j) - est_h(:,j);

    end
    
    for j = 1:max_iter
        % L2 cost in the SRVF domain
        for i = 1:N
           SRVF_L2cost(j) = SRVF_L2cost(j) + SquareL2norm(t, q_f(:,i) - MyGroupAction(t, est_q_g_plus_h(:,j),gam(:,i,j),3));
        end
        SRVF_L2cost(j) = SRVF_L2cost(j)/N;
        
        % L2 cost in the time domain
        for i = 1:N
            L2cost(j) = L2cost(j) + SquareL2norm(t, f(:,i) - MyGroupAction(t, est_g(:,j) + est_h(:,j), gam(:,i,j),1));
        end 
        L2cost(j) = L2cost(j)/N;
    end
    
    %%% compute error if synthetic data
    if true_data.use_data == 1
        g_error = L2norm(t, est_g(:,end) - true_data.g)/L2norm(t,true_data.g);   
        h_error = L2norm(t, est_h(:,end) - true_data.h)/L2norm(t,true_data.h);
        r_error = 0;
        for i = 1:N
            r_error = r_error + L2norm(t, gam(:,i,end) - true_data.r(:,i));
        end
        r_error = r_error/N;
        ortho = Check_Orthogonality(est_g(:,end),est_h(:,end), t);
    else
        ortho = Check_Orthogonality(est_g(:,end),est_h(:,end), linspace(0,1,true_data.T));
    end
    
    est_data.est_q_g_plus_h = est_q_g_plus_h;
    est_data.est_g_plus_h = est_g_plus_h;
    est_data.est_g = est_g;
    est_data.est_h = est_h;
    est_data.gam = gam;
    est_data.L2cost = L2cost;
    est_data.SRVF_L2cost = SRVF_L2cost;
    est_data.final_L2cost = L2cost(end);
    est_data.final_SRVF_L2cost = SRVF_L2cost(end);
    if true_data.use_data == 1
        est_data.g_error = g_error;
        est_data.h_error = h_error;
        est_data.r_error = r_error;
    end
    est_data.ortho = ortho;
end

