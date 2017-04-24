% This routine calculates the group action of gamma acting on f
% created at May, 10, 2014
% continued at June 4, 2015

function [ warp_f ] = MyGroupAction( t, f, gamma, choice )
% Input Argument:
%   t       of size Tx1
%   f       of size Tx1
%   gam     of size Tx1
%   choice  1 is f(r); 2 is f(r)*dr; 3 is f(r)*sqrt(dr) 
%
% Output Argument:
%   warp_f  of size Tx1, warped f

    if choice == 1
        %warp_f(:,1) = spline(t, f(:,1), gamma(:,1));
        warp_f(:,1) = interp1(t, f(:,1), gamma(:,1));
        
    elseif choice == 2
        [T,~] = size(t);  
        deri_f = gradient(gamma,1/(T-1));  
        %warp_f = spline(t,f,gamma).*(deri_f);
        warp_f = interp1(t,f,gamma).*(deri_f);
        
    elseif choice == 3
        [T,~] = size(t);  
        deri_f = gradient(gamma,1/(T-1));  
        %warp_f = spline(t,f,gamma).*sqrt(deri_f);
        warp_f = interp1(t,f,gamma).*sqrt(deri_f);
        
    else
        disp('wrong index in MyGroupAction');
    end
    
    if ( isnan(warp_f(1,1)) || isnan(warp_f(end,1)))
        disp('wrong in subroutine MyGroupAction');        
    end
end

