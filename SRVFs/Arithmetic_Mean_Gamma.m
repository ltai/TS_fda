% this is a program for computing arithmetic mean of given gamma
% created on Sep. 1, 2015;

function [ AM_gam ] = Arithmetic_Mean_Gamma( gam )

    [T,N] = size(gam);
    AM_gam = zeros(T,1);
    
    for i = 1:N
        AM_gam = AM_gam + gam(:,i);
    end
    AM_gam = AM_gam/N;

end

