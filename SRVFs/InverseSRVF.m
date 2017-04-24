% this is a routine for computing the inverse SRVF
% Warning: to gaurantee the unique transformation, require a initial point
% 
% created at Nov. 4, 2015; successful

function [ f ] = InverseSRVF( t,q,initial_pt )

    [T,~] = size(q);
    f = zeros(T,1);
    binsize = mean(diff(t));

    for j = 1:T
        if j == 1
            f(j) = 0;
        else
            temp_sum = 0;
            for k = 2:j
                temp_sum = temp_sum + ( q(k-1)*abs(q(k-1)) + q(k)*abs(q(k)) )*binsize/(2);
            end
            f(j) = temp_sum;
        end
    end

    f = f + initial_pt;
    
    
end

