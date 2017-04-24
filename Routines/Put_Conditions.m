% this program is for putting L1,L2 conditions on the input function
% so both g and h can be applied on the same way
% reduce the work of generate_g,h, cond_g,h

% created on July, 23, 2015

function [ fcn_out ] = Put_Conditions( fcn_in, para, idx_gh)
% idx_Norm is what type of L1,L2 conditions       
% idx_gh is for g or h; 1 is g, 2 is h 
        t = linspace(0,1,para.T)';        

        if idx_gh == 1
            fcn_L1norm = para.g_L1norm;  
            fcn_L2norm = para.g_L2norm;
            idx_Norm = para.g_cond;
              
        elseif idx_gh == 2
            fcn_L1norm = para.h_L1norm;
            fcn_L2norm = para.h_L2norm;
            idx_Norm = para.h_cond;
            
        end

        if idx_Norm == 0
            fcn_out(:,1) = fcn_in(:,1);
            
        elseif idx_Norm == 1
            fcn_out(:,1) = fcn_in(:,1) - trapz( t,fcn_in(:,1) ) + fcn_L1norm;
            
        elseif idx_Norm == 2
            fcn_out(:,1) = fcn_L2norm*fcn_in(:,1)/L2norm(t,fcn_in(:,1));
            
        elseif idx_Norm == 3
            %fcn_out(:,1) = fcn_in(:,1) + t*eps - trapz( t,fcn_in(:,1) ) + fcn_L1norm;
            fcn_out(:,1) = fcn_in(:,1) - trapz( t,fcn_in(:,1) ) + fcn_L1norm;
            fcn_out(:,1) = fcn_L2norm*fcn_out(:,1)/L2norm(t,fcn_out(:,1));
            
        elseif idx_Norm == 4
            fcn_out(:,1) = fcn_in(:,1) - trapz(t,fcn_in(:,1)) + 1;
            
        elseif idx_Norm == 100
            temp = fcn_in(:,1) - trapz(t,fcn_in(:,1));
            fcn_out(:,1) = fcn_L2norm*temp/L2norm(t,temp);

        else 
            fprintf('wrong in Put_Conditions');
        end
    
end

