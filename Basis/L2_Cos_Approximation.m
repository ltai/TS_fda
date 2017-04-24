% this program is doing L2 approximation between interval [0,1]
% using cosine basis.
% 
% Default: At most 100 basis. cos(n*pi*x), n = 0,1,2...
%                                      but in the program, index from 1
% created on July, 16, 2015; successful

function [ approx_f, Cos_Coef ] = L2_Cos_Approximation( f, start_idx, end_idx)
    
    NBasis = 100;
    
    if end_idx > NBasis
        warning('Exceed default max number of basis, L2_Cos_Approximation\n');
    end
    
    Cos_Coef = zeros(NBasis,1);
    
    [T,~] = size(f);
    t = linspace(0,1,T)';
    
    % form coefficient
    for k = start_idx:end_idx
        if k == 1
            Cos_Coef(k) = L2_InnerProduct(f,0*t+1,t);
        else
            Cos_Coef(k) = 2*L2_InnerProduct(f, cos((k-1)*pi*t),t);
        end
    end
    
    % form summation
    approx_f = zeros(T,1);
    
    for k = start_idx:end_idx
        approx_f = approx_f + Cos_Coef(k)*cos((k-1)*pi*t);
    end

end

