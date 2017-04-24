% this program is doing L2 approximation between interval [0,1]
% using cosine basis.
% 
% Default: At most 100 basis. sin(n*pi*x), n = 1,2,3...
%
% created on July, 16, 2015; successful
% continued on Feb. 16, 2016; double check

function [ approx_f, Sin_Coef ] = L2_Sin_Approximation( f, start_idx, end_idx)
    
    NBasis = 100; 
    if end_idx >= NBasis
        warning('Exceed default max number of basis, L2_Sin_Approximation\n');
    end
    Sin_Coef = zeros(NBasis,1);
    
    [T,~] = size(f);
    t = linspace(0,1,T)';
    
    % form coefficient
    for k = start_idx:end_idx
            Sin_Coef(k) = 2*L2_InnerProduct(f, sin(k*pi*t),t);
    end
    
    % form summation
    approx_f = zeros(T,1);
    
    for k = start_idx:end_idx
        approx_f = approx_f + Sin_Coef(k)*sin(k*pi*t);
    end

end

