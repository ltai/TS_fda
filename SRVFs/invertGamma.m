% don't understand this function in mathematics!!
% spline can be replaced by interp1 ?
% created at April, 29, 2014
% modified at May, 10, 2014: replace spline by interp1
% modified at May, 22, 2014: find the NaN problem in interp1
% condinued on Nov. 1, 2016; NaN problem still to be solved
%                soln: need odd T, otherwise return error
%                      also don't put the rescale to [0,1]

function gamI = invertGamma(gam)
% Input Argument
%   gam  of size (T x 1)
% Output Argument
%   gamI of size (T x 1), inverse of gamma

    version = 1;

    if version == 1     %%%% May,22,2014

        [T,~] = size(gam);
        t = linspace(0,1,T)';

        % gamI = spline(gam,x,x);  % spline will have problem at the boundary!
        % gamI = gamI/gamI(T);

        gamI = interp1(gam,t,t);
        %gamI = (gamI-min(gamI))/(max(gamI)-min(gamI)); 
        if (isnan(gamI(T,1)))
            gamI(T,1) = 1;
        end
        
        % ironically, interp1 is better than spline; Dec. 7, 2014 
    
    elseif version == 2 %%%%% May, 25, 2014, detect NaN and fix it
    
       [T,~] = size(gam);
       t = (0:T-1)'/(T-1);
       gamI = interp1(gam,t,t);
       gamI = gamI/gamI(T);
        if (isnan(gamI(1,1)))
            gamI(1,1) = 0;
        end
        
    elseif version == 3 %%%% Jan. 19, 2016, % same as version 1, but
                        % try to get rid of Nah
          [T,~] = size(gam);
          t = linspace(0,1,T)';
                        
          gamI = interp1(gam,t,t);
   
          for i = 1:T
              if (isnan(gamI(1,1)))
              end
          end
          gamI = (gamI-min(gamI))/ (max(gamI)-min(gamI));
                                         
    end    
end
