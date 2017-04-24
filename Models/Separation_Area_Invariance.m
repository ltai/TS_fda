% created on Jan. 7, 2016, need to revise
% continued on Jan. 13, 2016
% continued on Jan. 19, 2016; get stuck by the sqrt of (f_i-h),
%                             here require positive value!! finish code
%                             first
% comment: very careful about sqrt(f-h) and sqrt(g), that destroyed results
%

function [est_data] = Separation_Area_Invariance(para, true_data)

    T = true_data.T;
    N = para.N;
    f = true_data.f;
    max_iter = para.max_iter;
    t = true_data.t;
    
    h_start_idx = para.h_start_idx;
    h_end_idx = para.h_end_idx;
    
    est_g = zeros(T,max_iter);
    est_h = zeros(T,max_iter);
    gam = zeros(T,N,max_iter);
    gamI = zeros(T,N,max_iter);
    L2cost = zeros(max_iter,1);
    L2cost2 = zeros(max_iter,1);
       
    est_g(:,1) = mean(f,2);    
    for i = 1:N
        gam(:,i,1)  = t;
    end
    
    %%%% solve h and hereafter, h is known
    j = 2;
    for i = 1:N
        est_h(:,j) = est_h(:,j) + L2_Shift_Legendre_Approximation(f(:,i), h_start_idx, h_end_idx);
    end
    est_h(:,j) = est_h(:,j)/N;
    
    for j = 3:max_iter
        est_h(:,j) = est_h(:,2);
    end
    
    for j = 2:max_iter            
                
        %%%% solve gam
        lambda = 0;
        for i = 1:N
            gam(:,i,j) = MyDynamicProgramming( est_g(:,j-1)', sqrt(f(:,i) - est_h(:,j))' ,lambda,3)';
        end
        gam(:,:,j) = Gamma_with_inverse_kmean_id(gam(:,:,j));
        
        %%%% solve q_g and then g
        for i = 1:N
            gamI(:,i,j) = invertGamma(gam(:,i,j));  
            est_g(:,j) = est_g(:,j) + MyGroupAction(t, sqrt(f(:,i) - est_h(:,j)), gamI(:,i,j),3); 
        end
        est_g(:,j) = est_g(:,j)/N;   
        est_g(:,j) = est_g(:,j).^2;  
        
        est_g(:,j) = est_g(:,j) - L2_Shift_Legendre_Approximation(est_g(:,j), h_start_idx, h_end_idx);
        
    end
    
    for j = 1:max_iter        
        %%% L2 cost of l.h.s minus r.h.s
        for i = 1:N
            L2cost(j) = L2cost(j) + SquareL2norm(t, f(:,i) - MyGroupAction(t, est_g(:,j), gam(:,i,j),2) + est_h(:,j) );
        end 
        L2cost(j) = L2cost(j)/N;
        
        %%% L2 cost of the cost function
        for i = 1:N
            L2cost2(j) = L2cost2(j) + SquareL2norm(t, sqrt(f(:,i) - est_h(:,j)) ...
                - MyGroupAction(t,sqrt(est_g(:,j)), gam(:,i,j), 3) );
        end
        L2cost2(j) = L2cost2(j)/N;
    end

    est_data.est_g = est_g;
    est_data.est_h = est_h;
    est_data.gam = gam;
    est_data.L2cost = L2cost;
    est_data.L2cost2 = L2cost2;
end

