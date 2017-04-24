% this is a routine for computing Confidence interval
% created on July,16,2016; done :)

function [ upper_fcn,lower_fcn ] = Confidence_Interval( true_est_fcn, bs_fcn )
    [T,B] = size(bs_fcn);
    avg_bs_fcn = mean(bs_fcn,2);
    se_fcn_band = zeros(T,1);
    for grid = 1:T
        for b = 1:B
            se_fcn_band(grid) = se_fcn_band(grid) + (bs_fcn(grid,b)-avg_bs_fcn(grid)).^2;
        end
        se_fcn_band(grid) = se_fcn_band(grid)/(B-1);
        se_fcn_band(grid) = sqrt(se_fcn_band(grid));
    end
    upper_fcn = true_est_fcn + 1.96*se_fcn_band;
    lower_fcn = true_est_fcn - 1.96*se_fcn_band;
end

