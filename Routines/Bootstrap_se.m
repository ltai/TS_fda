% this routine computes bootstrap standard error
% created on July, 16, 2016; finished;

function [ se_replicate ] = Bootstrap_se( replicate )
    [~,B] = size(replicate);
    avg_replicate = mean(replicate);
    se_replicate = 0;
    for b = 1:B
        se_replicate = se_replicate + (replicate(b) - avg_replicate).^2;
    end
    se_replicate = se_replicate/(B-1);
    se_replicate = sqrt(se_replicate);   
end

