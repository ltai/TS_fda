% this is a routine for computing high order accuracy
% of the c.d.f. of standard normal distribution
% created on July,21,2016; successful;
% however, I found out this routine is no longer needed...

function [ cumulative_prob ] = My_cdf( given_value,mu,sigma )
    if given_value < -15 || given_value > 15
        printf('too exmreme value in My_cdf'); reture(0);
    end
    
    if nargin == 1
        mu = 0; sigma = 1;
    end
    T = 10^8;
    f = zeros(1,T);
    t = linspace(-100,given_value,T);
    cumulative_prob = 0;
    base = (given_value+100)/T;
    
    for i = 1:T
        f(i) = (1/sqrt(2*pi))*exp(-0.5*t(i)^2);
    end
    
    for i = 1:T-1
        cumulative_prob = cumulative_prob + 0.5*(f(i) + f(i+1))*base;
    end
end

