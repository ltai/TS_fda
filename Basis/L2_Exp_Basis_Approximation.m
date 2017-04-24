% this is a program for L2 of 1,e^a*t basis approximation,
% created on Oct. 15, 2015

function [ approx_f, approx_coef ] = L2_Exp_Basis_Approximation( f, start_idx, end_idx )

    Order = 10;  % 20 is best, otherwise loose numerical accuracy 
    if end_idx > Order
        warnnig('Exceed default max number of basis, L2_Exp_Basis_Approximation\n');
    end
    
    [T,~] = size(f);
    t = linspace(0,1,T)';
   
    %%% dense grid of sLegendre Poly.
    TT = 10^3+1;
    tt = linspace(0,1,TT)';
    % Generate_sLegendre_Poly(20,1001)
    
    temp = load('Ortho_Exp_Decay_Basis_1001','Ortho_Exp_Basis');
    Ortho_Exp_Basis = temp.Ortho_Exp_Basis;
    
    % L2 approximate f
    spline_f = spline(t,f,tt);
    
    approx_coef = zeros(Order,1);    % get approximation coefficient
    for k = start_idx:end_idx
        approx_coef(k) = (2*(k-1)+1)*L2_InnerProduct(spline_f,Ortho_Exp_Basis(:,k),tt);
    end
    
    approx_sf = zeros(TT,1);         % get approximation function
    for k = start_idx:end_idx
        approx_sf = approx_sf + approx_coef(k)*Ortho_Exp_Basis(:,k);
    end
    
    approx_f = spline(tt,approx_sf,t); % back to original grid numbers

end