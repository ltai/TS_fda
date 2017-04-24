% creatd on Jan.7 ,2016, 
% finished on Jan. 8, 2016, successful! Multiple iterations
% continued on Jan. 12, 2016, add L2 cost in the time domain

function [est_data] = Alignment_Model(para, true_data)

    T = true_data.T;
    N = true_data.N;
    f = true_data.f;
    max_iter = para.max_iter;
    t = true_data.t;
    
    est_q_g = zeros(T,max_iter);
    est_g = zeros(T,max_iter);
    gam = zeros(T,N,max_iter);
    gamI = zeros(T,N,max_iter);
    q_f = zeros(T,N);
    L2cost = zeros(max_iter,1);
    SRVF_L2cost = zeros(max_iter,1);
    
    %[fn,gam,qn,mq,gamI] = mainWarpingWrapper(t,f,0);
    
    for i = 1:N
        q_f(:,i) = SRVF(t,f(:,i));
    end
    est_q_g(:,1) = mean(q_f,2);
    
    for i = 1:N
        gam(:,i,1)  = t;
    end
    
    for j = 2:max_iter        
        %%%%% update gam
        lambda = 0;
        for i = 1:N
            gam(:,i,j) = MyDynamicProgramming( est_q_g(:,j-1)', q_f(:,i)',lambda,3)';
        end
        gam(:,:,j) = Gamma_with_inverse_kmean_id(gam(:,:,j));
        
        %%%%% update q_g        
        for i = 1:N
            gamI(:,i,j) = invertGamma(gam(:,i,j));     
            est_q_g(:,j) = est_q_g(:,j) + MyGroupAction(t, q_f(:,i), gamI(:,i,j),3);                
        end
        est_q_g(:,j) = est_q_g(:,j)/N;    
        
        est_g(:,j) = InverseSRVF(t,est_q_g(:,j),mean(f(1,:)));     
    end
    
    for j = 1:max_iter
        % L2 cost in the SRVF domain
        for i = 1:N
           SRVF_L2cost(j) = SRVF_L2cost(j) + SquareL2norm(t, q_f(:,i) - MyGroupAction(t, est_q_g(:,j),gam(:,i,j),3));
        end
        SRVF_L2cost(j) = SRVF_L2cost(j)/N;
        
        % L2 cost in the time domain
        for i = 1:N
            L2cost(j) = L2cost(j) + SquareL2norm(t, f(:,i) - MyGroupAction(t, est_g(:,j), gam(:,i,j),1));
        end 
        L2cost(j) = L2cost(j)/N;
    end
    
    if true_data.use_data == 1
        g_error = L2norm(t, est_g(:,end) - true_data.g)/L2norm(t,true_data.g);   
        r_error = 0;
        for i = 1:N
            r_error = r_error + L2norm(t, gam(:,i,end) - true_data.r(:,i));
        end
        r_error = r_error/N;
    end
    
    est_data.est_q_g = est_q_g;
    est_data.est_g = est_g;
    est_data.gam = gam;
    est_data.L2cost = L2cost;
    est_data.SRVF_L2cost = SRVF_L2cost;
    if true_data.use_data == 1
        est_data.g_error = g_error;
        est_data.r_error = r_error;
    end
end
