% not yet finished
%
%  Default: At most 100 basis. sin(2n*pi*x) + cos(2*n*pi*x), n = 0,1,2...
%  At most 200 Basis 
%
% created on July, 23, 2015;
% continued on July, 29, 2015; finished

function [ approx_f, MFS_Coef ] = L2_MFS_Approximation( f, start_idx, end_idx )
    [T,~] = size(f);
    t = linspace(0,1,T)';
    
    NBasis = 200;
    
    if end_idx >= NBasis/2
        warning('Exceed default max number of basis, L2_MFS_Approximation\n');
    end
    half_NBasis = NBasis/2;
    Cos_Coef = zeros(half_NBasis,1);
    Sin_Coef = zeros(half_NBasis,1);

    % get correct coefficient
    for k = start_idx:end_idx
        if k == 1
            Cos_Coef(k) = L2_InnerProduct(f, cos(2*pi*(k-1)*t), t);
            Sin_Coef(k) = 0;
        else
            Cos_Coef(k) = 2*L2_InnerProduct(f, cos(2*pi*(k-1)*t), t);
            Sin_Coef(k) = 2*L2_InnerProduct(f, sin(2*pi*(k-1)*t), t);
        end
    end

    MFS_Coef = [Cos_Coef'; Sin_Coef']';
    
    % use coefficient to form f, L2 approximation
    approx_f = zeros(T,1);
    for k = start_idx:end_idx
        approx_f = approx_f + Cos_Coef(k).*cos(2*pi*(k-1)*t)...
                + Sin_Coef(k).*sin(2*pi*(k-1)*t);
    end
   
end

