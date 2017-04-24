% approximated L2 in [0,1] by using shift Legendre polynomial
% created on Aug. 5, 2015
% revised original one into approximated spline f, dense points to
%                  reduce numerical error a lot!
% continued on Aug. 6, 2015; success; speed up.
% continued on Aug. 13, 2015; 
%                (idx1) spline slow donw but high L2 approximation accuracy
%                (idx2) Simpson is fast but low L2 approximation accuracy    

function [ approx_f, approx_coef ] = L2_Shift_Legendre_Approximation( f, start_idx, end_idx )
    
    Order = 20;  % 20 is best, otherwise loose numerical accuracy 
    if end_idx > Order
        warnnig('Exceed default max number of basis, L2_sLegendre_Approximation\n');
    end
    
    [T,~] = size(f);
    t = linspace(0,1,T)';
        
    idx = 1; % 1 is dense grid; 2 is coarse grid
    
    if idx == 1               %%% dense grid of sLegendre Poly.
        TT = 10^3+1;
        tt = linspace(0,1,TT)';
        % Generate_sLegendre_Poly(20,1001)

        temp = load('Legendre_Poly_1001','Legendre_Poly');
        Legendre_Poly = temp.Legendre_Poly;

        % L2 approximate f
        spline_f = spline(t,f,tt);    
        
        approx_coef = zeros(Order,1);    % get approximation coefficient
        for k = start_idx:end_idx
            approx_coef(k) = (2*(k-1)+1)*L2_InnerProduct(spline_f,Legendre_Poly(:,k),tt);
        end
        
        approx_sf = zeros(TT,1);         % get approximation function
        for k = start_idx:end_idx
            approx_sf = approx_sf + approx_coef(k)*Legendre_Poly(:,k);
        end

        approx_f = spline(tt,approx_sf,t); % back to original grid numbers
    
    elseif idx == 2    %%% coarse grid of sLegendre Poly.      
        temp = load('Legendre_Poly_101','Legendre_Poly');
        Legendre_Poly = temp.Legendre_Poly;
    
        approx_coef = zeros(Order,1);
        
        for k = start_idx:end_idx
            approx_coef(k) = (2*(k-1)+1)*L2_InnerProduct(f,Legendre_Poly(:,k),t);
        end
        
        approx_f = zeros(T,1);         % get approximation
        for k = start_idx:end_idx
            approx_f = approx_f + approx_coef(k)*Legendre_Poly(:,k);
        end
        
    end
 end

