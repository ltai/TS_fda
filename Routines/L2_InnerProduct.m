% This is a routine computing L2 inner product using trapezoidal rules
% create at Dec. 7, 2014; working
% continued on Aug. 5, 2015; idx = 2 is very accurate since involving 10^5
%                            points; however, very slow since need spline

function [ L2_IP ] = L2_InnerProduct( f, g, t )
% Input:
%    f,g are two input functions, of size Tx1 
%    t   is time grid, of siez Tx1      
% Output:
%    L2_IP is the L2 inner product
    idx = 3;
    
    if idx == 1                       % old version
        [T,~] = size(t);
        product = f.*g;
        L2_IP = 0;

        for i = 1:T-1
            L2_IP = L2_IP + ( product(i) + product(i+1)) * ( t(i+1)-t(i) )/2;
        end
        
    elseif idx == 2           % new version, add more points and with spline
        TT = 10^4;
        tt = linspace(0,1,TT)';
        
        ff = spline(t,f,tt);
        gg = spline(t,g,tt);
        
        L2_IP = trapz(tt,ff.*gg);
    
    elseif idx == 3          % trapezoidal rule to evaluate quadrature
        
        a = 0; b = 1;
        L2_IP = Simpsons_Rule( f.*g ,a,b );
        
    end
end

