% This is the Algorithm 1 of paper Srivastava, 2011,
%                         Karcher mean of points in Gamma
% created by April, 29, 2014; this is the original version
% testing on Sep. 6, 2014
% continued on Sep. 1, 2015; data structure need to be revised.

function [KM_gam] = Karcher_Mean_Gamma(gam)
%
% Input Argument
%   gam of size (T x N); N is number of gammas; T is time steps
%
% Output Argument
%   gam of size (Tx1);
  
    gam = gam';          % still didn't change the data structure
    [N,T] = size(gam);
    
    dT = 1/(T-1);
    psi = zeros(N,T-1);
    
    for i=1:N
        psi(i,:) = sqrt( diff(gam(i,:))/dT ); % lose one grid point
        %psi(I,:) = sqrt( gradient(gam(i,:),dT) );
    end

    % Find direction
    mu = psi(1,:);
    vec = zeros(N,T-1);
    t = 1;
    clear vec;
    
    for iter = 1:50 % This step is strange ! why choose step # to be 10 only
        for i = 1:N
            %v = psi(i,:) - mu;
            len = acos(sum(mu.*psi(i,:))*dT);
            
            if len > 10^-8            
                vec(i,:) = (len/sin(len))*(psi(i,:) - cos(len)*mu);
            else
                vec(i,:) = zeros(1,T-1);
            end        
        end
        vm = mean(vec);
        
        lvm(iter) = sqrt(sum(vm.*vm)*dT);
        if lvm  > 10^-8
            mu = cos(t*lvm(iter))*mu + (sin(t*lvm(iter))/lvm(iter))*vm;
        end
    end

    gam_mu = [0 cumsum(mu.*mu)]/T;
    gam_mu = (gam_mu-min(gam_mu))/(max(gam_mu)-min(gam_mu));
    
    KM_gam = gam_mu'; % want ourput a column
   
end
