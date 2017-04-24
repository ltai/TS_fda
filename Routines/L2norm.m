% This routine gives the L^2 norm of a function on interval [0,1]
% by using quadrature rule

% created at May, 11, 2014
% modified on Sep. 8, 2014; change number of grid points depends
%                           on the input t;
% comment: using linear interpolation won't help in trapezondal rule
%          instead, try spline!!
% continued on Aug. 12, 2015; revised to spline interpolation

function [ L2norm ] = L2norm( t,f )
% Input Argument
%   t        of size Tx1
%   f        of size Tx1
% Output Argument
%   L2norm   scalar

    idx = 3;
    
    if idx == 1 % odd version; can be deleted
        
        [T,~] = size(t);
        Npts = 5*T; % number of timesteps used in the quadrature
        L2norm = 0;
        
        tt = linspace(0,1,Npts)';
        %ff = interp1(t,f,tt);
        ff = spline(t,f,tt);
        
        if( isnan(ff(1,1)) || isnan(ff(end,1)))
            warning('wrong in subroutine L2_norm\n');
        end
        
        %ff = spline(t,f,tt); % didn't see difference using interp1 and spline
        ff = ff.^2;
        
        for i = 1:Npts-1
            L2norm = L2norm + ( ff(i) + ff(i+1) )*( tt(i+1) - tt(i) )/2;
        end
        L2norm = sqrt(L2norm);
    
    elseif idx == 2  % new version; can deleted
        TT = 10^4;
        tt = linspace(0,1,TT)';
        ff = spline(t,f,tt);
        
        SquareL2norm = trapz(tt,ff.*ff);
        L2norm = sqrt(SquareL2norm);
    
     elseif idx == 3 %%%%%   Simpson's Rule    %%%%%%
        a = 0; b = 1;
        quad = Simpsons_Rule( f.*f ,a,b );
        L2norm = sqrt(quad);
        
    else
        warning('wrong idx in L2nrom\n');
    end
    
end

