% This is a routine for computing gradient. In the Matlab comment, gradient
% the middle point has O(h^2) but the end points only has O(h)
% Now the derivative use O(h^4) for middle points (5-point stencle)
% and O(h^2) for 1,2, T-1,T points.

% created at June, 11, 2014

 function [dev_f] = MyGradientV2(f,t)
% Input
%       t   of size Tx1, the time grid
%       f   of size Tx1, function value, already in vector form
% Output
%       dev_f of size Tx1, the derivative of f
    %T = 100;
    %t = linspace(0,1,T)';
    % f = sin(2*pi*t);
    
    [T,~] = size(t);
    h = (t(end)-t(1))/(T-1);  % assume uniform mesh
    dev_f = zeros(T,1);
    
    if (T<5)
        disp('need more than 5 points when taking diff.');
        disp('wrong in MyGradient routine');
        return;
    end
    
    for i = 1:T 
        if (i == 1)
            %dev_f(i,1) = ( -3*f(i,1) + 4*f(i+1,1) -f(i+2,1) )/(2*h); % O(h^2)
            dev_f(i,1) = ( -21*f(i,1) +32*f(i+1,1) -12*f(i+2,1) +f(i+4,1) )/(12*h); % O(h^3)
            
        elseif (i == T)
            %dev_f(i,1) = ( -3*f(i,1) + 4*f(i-1,1) - f(i-2,1) )/(-2*h); % O(h^2)
            dev_f(i,1) = ( -21*f(i,1) +32*f(i-1) -12*f(i-2) +f(i-4) )/(-12*h); % O(h^3)
            
        elseif (i == 2 )
            %dev_f(i,1) = (f(i+1,1)-f(i-1,1))/(2*h);    % O(h^2)
            dev_f(i,1) = ( -2*f(i-1,1) -3*f(i,1) +6*f(i+1,1) -f(i+2,1) )/(6*h); % O(h^3)
            
        elseif (i == T-1)   
            %dev_f(i,1) = (f(i+1,1)-f(i-1,1))/(2*h);    % O(h^2)
            dev_f(i,1) = ( 2*f(i+1,1) +3*f(i,1) -6*f(i-1,1) +f(i-2,1) )/(6*h); % O()h^3
            
        else
            % O(h^4)
            dev_f(i,1) = ( f(i+1,1) - f(i-1,1) )/(2*h) ...
                       - ( f(i+2,1) - 2*(f(i+1,1)-f(i-1,1)) -f(i-2,1))/(12*h);
        end
        
%         if (dev_f(i,1)<10^-7)
%             disp('small magnitude in derivative!');
%         end
    end

    %df = inline('2*pi*cos(2*pi*t)','t');
    %dfcn = df(t);
    %nu_df = gradient(f,1/(T-1));
    %plot(t,dfcn,t,dev_f);
    
 end
