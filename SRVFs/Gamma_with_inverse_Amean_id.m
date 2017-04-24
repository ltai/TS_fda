% This routine is imposing inverse Arithmetic mean of inverse gamma equals
%                             identity
% created at Sep. 1, 2015

function [ imposed_gam,mu, muI] = Gamma_with_inverse_Amean_id( gam )
    
    [T,N] = size(gam);
    t = linspace(0,1,T)';
    imposed_gam = zeros(T,N);
    gamI = zeros(T,N);

    for i = 1:N
        gamI(:,i) = invertGamma(gam(:,i));
    end

    mu = Arithmetic_Mean_Gamma(gamI);
    muI = invertGamma(mu);
    
    for i = 1:N
        imposed_gam(:,i) = interp1(t,mu,gam(:,i));
    end
    
end

