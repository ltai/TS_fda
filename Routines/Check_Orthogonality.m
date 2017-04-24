% This is a routine computing L2 inner product using trapezoidal rules
% create at Dec. 7, 2014; working
% continued on Aug. 5, 2015; idx = 2 is very accurate since involving 10^5
%                            points; however, very slow since need spline
% continued on Aug. 20, 2015; replace quad to Simpson's Rule

function [ L2_IP ] = Check_Orthogonality( f, g, t )
% Input:
%    f,g are two input functions, of size Tx1 
%    t   is time grid, of siez Tx1      
% Output:
%    L2_IP is the L2 inner product

    idx = 2;
    
    if idx == 1               % directly using Simpson's Rule
        L2_IP = Simpsons_Rule(f.*g,0,1);
 
    elseif idx == 2           % dense spline with Simpson's Rule
        TT = 10^3+1;
        tt = linspace(0,1,TT)';
        
        ff = spline(t,f,tt);
        gg = spline(t,g,tt);
        
        L2_IP = Simpsons_Rule(ff.*gg,0,1);
        
    end
end

