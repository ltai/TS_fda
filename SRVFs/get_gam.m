% created at June, 11, 2014
% continud on June, 14, 2015
% continued on Nov. 1., 2016; done; successful :)

%%% choice 1
% generate gamma with r = t^a where a is a number not too far from 1
% time grid is assuming uniform from 0 to 1, total T points
% note the first point may be not defined, not yet shifted (ignore)
%
%%% choice 2
% this is a shited version of r(t) = t.^a 
% solves the undefined problem whtn t=0.
%
%%% choice 3
% idea from one of Tucker's code, but I adapt to my own version.
% he forget to put square if want a more general case
%
%%% choice 4
% Based on nonlinear transformation to a exponential function
% This gamma has the infinitely differentiability
% idea come from Srivastava, 2011, paper
%
%%% choice 5
% r = integrate of (cos)^2
%
%%% choice 6 (similar to 5)
% r = integrate of (cos+sin)^2

function [ r,rI ] = get_gam(T, N, r_choice)
% T       is time  
% choice  is the choice for generating gamma
% r       is gamma,
% rI      is inverse of gamma
    r = zeros(T,N);
    rI = zeros(T,N);
    t = linspace(0,1,T)';

    if r_choice == 1 
        p = 1;      % power, [0.5, 2] is suitable for practice
        r = t.^p;
        rI = t.^(1/p);
    
    elseif r_choice == 2
        shift = 1; % should away from 0, positive number

        if ( mod(N,2)==1 ) % odd N
            p1 = linspace(0.05, 0.8, (N-1)/2);   % [0.5, 2] is suitable
            p2 = linspace(1.1, 2, (N-1)/2);
            p = [p1 1 p2];
            
        elseif ( mod(N,2)==0 ) % even N
            p1 = linspace(0.05, 0.8, N/2);   % [0.5, 2] is suitable
            p2 = linspace(1.1, 2, N/2);
            p = [p1 p2];
        end
        
        for i = 1:N
            r(:,i) = ( (t'+shift).^p(i) - shift^p(i) )/( (1+shift)^p(i) - shift^p(i) );
            rI(:,i) = nthroot( ((1+shift)^p(i)-shift^p(i)).*t' + shift^p(i) , p(i)) - shift;
        end       
        
    elseif r_choice == 3     
        f = zeros(T,N);
        temp = linspace(-0.5,0.5,N);
        for i = 1:N
            % f(:,i) = sin(temp(i)*pi*t)+0.1;
             f(:,i) = cos(temp(i)*pi*t);
            % f(:,i) = -sin(temp(i)*pi*t) + 2*cos(temp(6-i)*pi*t)+0.1;
            %f(:,i) = f(:,i).^2+0.1;
            r(:,i) = cumtrapz(t,f(:,i))/trapz(t,f(:,i));
        end
        
    elseif r_choice == 4    %%%%%%% I use this one most of the time
        deri_r = zeros(T,N);
        deri_rI = zeros(T,N);
        %a = linspace(-1,8,N); % parameters (-3,3)
        a = linspace(-3,3,N); % (-2,2)
        %a = [-2 -1 2 4 6];
        %a = [-2 -1 1 3 5];       
        for i = 1:N
            if a(i) == 0
                r(:,i) = linspace(0,1,T)';
                rI(:,i) = linspace(0,1,T)';
                deri_r(:,i) = ones(T,1);
                deri_rI(:,i) = ones(T,1);
            else
                r(:,i) = ( exp(a(i)*t)-1 )/( exp(a(i))-1 );
                rI(:,i) = (1/a(i))*log( ( exp(a(i))-1 ).*t + 1 );
                deri_r(:,i) = ( a(i)*exp(a(i)*t) )/( exp(a(i))-1 );
                deri_rI(:,i) = ( exp(a(i))-1 )./( a(i)*( (exp(a(i))-1)*t+1 ) );
            end
        end
        
    elseif r_choice == 5       
        f = zeros(T,N);
        temp = linspace(-0.5,0.5,N);
        for i = 1:N
            f(:,i) = 3*cos(pi*t+temp(i));
            f(:,i) = f(:,i).^2+0.1;
            r(:,i) = cumtrapz(t,f(:,i))/trapz(t,f(:,i));
        end
        
    elseif  r_choice == 6
        f = zeros(T,N);
        temp = linspace(-2,2,N);
        for i = 1:N
            f(:,i) = 3*sin(temp(i)*2*pi*t)+2*cos(t*temp(i));
            f(:,i) = f(:,i).^2+1;
            r(:,i) = cumtrapz(t,f(:,i))/trapz(t,f(:,i));
        end   
    end    
end

