% created on July, 1,2016; successful;

function [ para ] = Set_Algorithm_Para(max_iter,update_h,h_start_idx,h_end_idx)
    
%%%%%%%%%%%%%%   set general parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if nargin == 0
    para.max_iter = 5;   % this is for Coordinate Descent  
    para.update_h = 30;
    para.h_start_idx = 1;   % start_idx of basis of h 
    para.h_end_idx = 4;     % end_idx of basis of h 
    
elseif nargin == 4
    para.max_iter = max_iter;
    para.update_h = update_h;
    para.h_start_idx = h_start_idx;
    para.h_end_idx = h_end_idx;
else
    fprintf('wrong number inputs of Set_Algorithm_Para'); return;
end
      
    para.g_start_idx = 5;   % start_idx of basis of g 
    para.g_end_idx = 20;    % end_idx of basis of g 
     
    %%%%%% choose the way of updating r,g, and h 
    para.update_r = 3;   % 0 is not update, true;
                         % 2 is update by group action f(r)*(dr); 
                         % 3 is update by group action f(r)*sqrt(dr);
    
    para.update_g = 40;  % 0 is not update;
                         % 1 is non-para average update; 
                         % 15 is para average, Cos Basis;
                         % 20 is para average, Sine Basis
                         % 25 is para average, Modified Fourier Basis
                         % 30 is para average, Shifted Legendre Basis
                         % 40 is non-para average with ortho. to para. h
                                                  
    %para.update_h = 30;  % 0 is not update
                         % 1 is non-para average update; 
                         % 10 is para average, Newton;
                         
                         % 15 is para average, Cos Basis;
                         % 20 is para average, Sine Basis
                         % 25 is para average, Modified Fourier Basis
                         % 30 is para average, Shifted Legendre Basis
                         
                         % 40 is non-para average with ortho. to para. g
            
%%%%%%%%%%%%%% no need to change below %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% no need to change below %%%%%%%%%%%%%%%%
    para.GA_choice = 3;     % group action choice; 2 is g(r)*dr; 
                            %                      3 is g(r)*sqrt(dr);
    para.est_r_cond = 1;    % 0 is no cond.;
                            % 1 is KM{r^-1 _i} = r_id;
                            % 2 is AM{r^-1 _i} = r_id;
    % no longer exist in the code                  
    para.est_g_cond = 0;    % 0 is no cond.; 
                            % 1 is int = 0; 
                            % 2 is L2-norm = 1; 
                            % 3 is int = 0 and L2-norm = 1; 
                            % 4 is int = 1; 
                            % 5 is from invariance % forget the purpose
                            % 100 is int(g) = 0 and L2norm(g) = c1
    % no longer exist in the code              
    para.est_h_cond = 0;    % 0 is no cond.; 
                            % 1 is int = 0; 
                            % 2 is L2-norm = 0.2; 
                            % 3 is int = 0 and L2-norm = 0.2; 
                            % 4 is L2-norm = 5; 
                            % 5 is int(abs) = 1;   
                            % 100 is int(h) = 0 and L2norm(h) = c2
                            % 101 is int(h) = c1 and L2norm(h) = c2; relate to f_i, g, and gam   
end

