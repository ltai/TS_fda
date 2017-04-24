% this routine gives the Square of L^2 norm of a function on interval [0,1]
% by using mid-point quadrature rule
%
% created at July, 20, 2014
% continued on Aug. 12, 2015; revised to spline interpolation,
%                             add Simpson's Rule

function [ SquareL2norm ] = SquareL2norm( t,f )
% Input Argument
%   t        of size Tx1
%   f        of size Tx1
% Output Argument
%   SquareL2norm   scalar

    idx = 3;
    
    if idx == 1 %%%%%%%%%%%% trapezondal rule  %%%%%%%%%%
        [T,~] = size(f);
        Npts = T; % number of timesteps used in the quadrature
        L2norm = 0;

        tt = linspace(0,1,Npts)';
        ff = interp1(t,f,tt);
        %ff = spline(t,f,tt);

        if( isnan(ff(1,1)) || isnan(ff(end,1)))
               disp('wrong in subroutine L2_norm');
        end

        %ff = spline(t,f,tt); % didn't see difference using interp1 and spline
        ff = ff.^2;

        for i=1:Npts-1
            L2norm = L2norm + ( ff(i) + ff(i+1) )*( tt(i+1) - tt(i) )/2;
        end
        SquareL2norm = L2norm;
    
    elseif idx == 2  %%%%%%  spline with trapezoidal     %%%%%%%%%
        TT = 10^4;
        tt = linspace(0,1,TT)';
        ff = spline(t,f,tt);
        
        SquareL2norm = trapz(tt,ff.*ff);
        
    elseif idx == 3 %%%%%   Simpson's Rule    %%%%%%
        a = 0; b = 1;
        SquareL2norm = Simpsons_Rule( f.*f ,a,b );

    else
        warning('wrong idx in SquareL2norm\n');        
    end
end

