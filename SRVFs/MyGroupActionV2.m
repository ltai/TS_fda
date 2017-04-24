% This routine computes group action:  (for)*sqrt(dot(r))
% Inproves previous routine "warp_f_by_gammaV3"
% The main difference is deri_rI should be replaced by 1/deri_r  

% comment: for boundary blow up, it seems not related to use spline or interp1

% created on Sep. 20, 2014; testing; O.K.
% continued on Nov. 2, 2014; theoretically, MyGradient has 6-th order
%                             accuracy, should have better result.
%                           However, MyGradient seems sensitive from 
%                           result gam = DynamicProgrammingV3.c
% continued on June, 3, 2015; This MyGroupActionV2 but has nothing to do
%                               with the other routine MyGroupAction.


function [ warp_f ] = MyGroupActionV2( t, f, gam, gamI, idx )
% Input Argument:
%   t       of size Tx1
%   f       of size Tx1
%   gam     of size Tx1
%   gamI    of size Tx1
%   idx     of size 1x1; 1 for (f,gam) and 2 for (f,gamI)
% Output Argument:
%   warp_f  of size Tx1, group action of f
    
    [T,~] = size(t);
    warp_f = zeros(T,1);
    
    if idx == 1
        %deri_gam = MyGradient(gam,t);
        deri_gam = gradient(gam,t);
        
        if deri_gam < 0 
            printf('wrong in GroupAction'); 
        end
        
        warp_f = spline(t,f,gam).*sqrt(deri_gam);
        %warp_f = interp1(t,f,gam).*sqrt(deri_gam);
        
    elseif idx == 2
        % Version 1
        % deri_gam = MyGradient(gam,t);
        % deri_gamI = 1./deri_gam;
        % deri_gamI = spline(gam,deri_gamI,t);
        
        % Version 2
        %deri_gamI = MyGradient(gamI,t);
        
        % Version 3
        %deri_gamI = gradient(gamI,t);
        
        % Version 4
         deri_gam = gradient(gam,t);
         deri_gamI = 1./deri_gam;
         deri_gamI = spline(gam,deri_gamI,t);
        
        if deri_gamI < 0 
            printf('wrong in GroupAction'); 
        end
        
        warp_f = spline(t,f,gamI).*sqrt(deri_gamI);
        %warp_f = interp1(t,f,gamI).*sqrt(deri_gamI);
    end
 
    if ( isnan(warp_f(1,1)) || isnan(warp_f(end,1)))
        disp('wrong: get NaN');        
    end
end

