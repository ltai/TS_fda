% created on Jan. 7, 2016
% finished on Jan. 8, 2016

function [est_data] = Simple_Estimation(para, true_data)
    
    T = true_data.T;
    N = true_data.N;
    f = true_data.f;
    max_iter = para.max_iter;
    t = true_data.t;
    
    est_g = zeros(T,max_iter);
    L2cost = zeros(max_iter,1);

    for j = 1:max_iter
        est_g(:,j) = mean(f,2);
    end

    for j = 1:max_iter
        temp_L2cost = 0;
        for i = 1:N
            temp_L2cost = temp_L2cost + SquareL2norm(t,f(:,i) - est_g(:,j));
        end
        L2cost(j) = temp_L2cost/N;
    end
    
    if true_data.use_data == 1
        g_error = L2norm(t, est_g(:,end) - true_data.g)/L2norm(t,true_data.g);   
    end
    
    est_data.est_g = est_g;
    est_data.L2cost = L2cost;
    if true_data.use_data == 1
        est_data.g_error = g_error;
    end
end

