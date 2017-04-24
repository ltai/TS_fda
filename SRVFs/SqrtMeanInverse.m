% Modified at May, 10, 2014; want to chage the data structure of gamma
%                            but indeed, not yet

function gamI = SqrtMeanInverse(gam)
% This is the Algorithm 1, Karcher mean of points in Gamma
%                          and find it's inverse
% Input Argument
%   gam of size T x N; T is # of time steps; N is # of gammas
% Output Argument
%   gamI of size T x 1
%
    gam = gam';
    [N,T] = size(gam);
    dT = 1/(T-1);
    psi = zeros(N,T-1);
    for i=1:N
        psi(i,:) = sqrt(diff(gam(i,:))/dT);
    end

    %Find direction
    mu = psi(1,:);
    t = 1;
    clear vec;
    for iter = 1:10 % This step is strange ! why choose step # to be 5 only
        for i=1:N
            v = psi(i,:) - mu;
            len = acos(sum(mu.*psi(i,:))*dT);
            if len > 0.000001            
                vec(i,:) = (len/sin(len))*(psi(i,:) - cos(len)*mu);
            else
                vec(i,:) = zeros(1,T-1);
            end        
        end
        vm = mean(vec);
        lvm(iter) = sqrt(sum(vm.*vm)*dT);
        if lvm  > 0.000001
            mu = cos(t*lvm(iter))*mu + (sin(t*lvm(iter))/lvm(iter))*vm;
        end
    end

    gam_mu = [0 cumsum(mu.*mu)]/T;
    gam_mu = (gam_mu-min(gam_mu))/(max(gam_mu)-min(gam_mu));
    gamI = invertGamma(gam_mu');
    
end
