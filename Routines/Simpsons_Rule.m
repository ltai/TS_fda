% this is a routine for computing numerical quadrature
% by using Simpson's Rule
%
% created on Aug. 13, 2015; successful

function [ quad ] = Simpsons_Rule( temp_f,a,b )
    
    [temp_T,~] = size(temp_f);
    
    if mod(temp_T,2) == 0
        warning('need odd grid points in SquareL2norm, re_spline');
        T = temp_T + 1;
        t = linspace(a,b,temp_T);
        tt = linspace(a,b,T);
        f = spline(t,temp_f,tt);
        
    else
        f = temp_f;
        T = temp_T;
    end
    quad_temp = 0;
    
    for i = 1:T
        
        if i == 1 || i == T
            quad_temp = quad_temp + f(i);
            
        elseif mod(i,2) == 0
            quad_temp = quad_temp + 4*f(i);
            
        elseif mod(i,2) == 1
            quad_temp = quad_temp + 2*f(i);
            
        else
            warning('wrong in Simpsons_Rule\n');
        end
    end
    quad = quad_temp * (b-a)/(3*(T-1));

end

