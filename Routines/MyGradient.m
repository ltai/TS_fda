% This is a routine for computing gradient.
% Improve the original MyGradient which only have O(h^3) on the 
% boundary and O(h^4) in the middle
% All are O(h^-6)

% created at Sep, 8, 2014; working :)
% modified on Sep. 20, 2014; go to O(h^-8) ! Working.

 function [dev_f] = MyGradient(f,t)
% Input
%       t   of size Tx1, the time grid
%       f   of size Tx1, function value, already in vector form
% Output
%       dev_f of size Tx1, the derivative of f

    
    [T,~] = size(t);
    h = (t(end)-t(1))/(T-1);  % assume uniform mesh
    dev_f = zeros(T,1);
    order = 6;
    
    if (T<7)
        disp('need more than 7 points when taking diff.');
        disp('wrong in MyGradient routine');
        return;
    end
    
    if order == 8
        for i = 1:T
            if i<=4
                dev_f(i,1) = (-761/280)*f(i,1) + (8)*f(i+1,1) + (-14)*f(i+2,1) + (56/3)*f(i+3,1) +...
                    (-35/2)*f(i+4,1) + (56/5)*f(i+5,1) + (-14/3)*f(i+6,1) + (8/7)*f(i+7,1) + (-1/8)*f(i+8,1);
                dev_f(i,1) = dev_f(i,1)/h;
                
            elseif T-3 <= i
                dev_f(i,1) = (1/8)*f(i,1) + (-8/7)*f(i-1,1) + (14/3)*f(i-2,1) + (-56/5)*f(i-3,1) +...
                    (35/2)*f(i-4,1) + (-56/3)*f(i-5,1) + (14)*f(i-6,1) + (-8)*f(i-7,1) + (761/280)*f(i-8,1);
                dev_f(i,1) = dev_f(i,1)/h;
                
            else
                dev_f(i,1) = (0)*f(i,1) + (4/5)*f(i+1,1) + (-1/5)*f(i+2,1) + (4/105)*f(i+3,1) +...
                    (-1/280)*f(i+4,1) + (-4/5)*f(i-1,1) + (1/5)*f(i-2,1) + (-4/105)*f(i-3,1) + (1/280)*f(i-4,1);
                dev_f(i,1) = dev_f(i,1)/h;
                
            end
        end
        
    elseif order == 6      
        for i = 1:T
            if (i <= 3)
                
                dev_f(i,1) = (-49/20)*f(i,1) + (6)*f(i+1,1) + (-15/2)*f(i+2,1) +...
                    (20/3)*f(i+3,1) + (-15/4)*f(i+4,1) + (6/5)*f(i+5,1) + (-1/6)*f(i+6,1) ;
                
                dev_f(i,1) = dev_f(i,1)/h;
                
            elseif ( T-2 <= i )
                dev_f(i,1) = (49/20)*f(i,1) + (-6)*f(i-1,1) + (15/2)*f(i-2,1) +...
                    (-20/3)*f(i-3,1) + (15/4)*f(i-4,1) + (-6/5)*f(i-5,1) + (1/6)*f(i-6,1) ;
                
                dev_f(i,1) = dev_f(i,1)/h;
            else
                
                dev_f(i,1) = (0)*f(i,1) + (3/4)*f(i+1,1) + (-3/20)*f(i+2,1) +...
                    (1/60)*f(i+3,1) + (-3/4)*f(i-1,1) + (3/20)*f(i-2,1) + (-1/60)*f(i-3,1);
                
                dev_f(i,1) = dev_f(i,1)/h;
            end
        end
    end  
        
    
      
 end
