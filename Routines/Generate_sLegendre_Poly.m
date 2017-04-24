% this is a program for generating shifted legendre polynomial
% reduce the time consuming of the orther routine:
%      L2_Shift_Legendre_ApproximationV2 
%
% created on Aug. 5, 2015; successful

%function  Generate_sLegendre_Poly(Order,TT)

    Order = 20;
    %TT = 10^4+1;
    TT = 101;
    tt = linspace(0,1,TT)';

    sLegendre_Coef = Shift_Legendre_Coef( Order );
    Legendre_Poly = zeros(TT,Order);

    % form normalized shift Legendre poly

    for n = 1:Order
        for k = 1:Order
            Legendre_Poly(:,n) = Legendre_Poly(:,n) + sLegendre_Coef(k,n)*tt.^(k-1);
        end
    end

    save('Legendre_Poly.mat','Legendre_Poly');
%end

