% this is generate coefficient of N-th order Legender polynomial 
% in [0,1],
% 
% created on July 16, 2015; successful

function [ Legendre_Coef ] = Shift_Legendre_Coef( Order )

    %Order = 50;   % Max Order is 85, otherwise out of precision

    Legendre_Coef = zeros(Order,Order);

    for n = 1:Order
        for k = 1:n

            comb1 = factorial(n-1)/( factorial(k-1)*factorial(n-k) );
            comb2 = factorial(n-1+k-1)/( factorial(k-1)*factorial(n-1) );

            Legendre_Coef(k,n) = Legendre_Coef(k,n) + (-1)^(n-1) * comb1 * comb2 * (-1)^(k-1);

        end
    end

end

