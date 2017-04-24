% this is a program for computing loss function, or
% cost function
%
% the "idx" represent using no regularization 0,1,10,15,20,25,30 => 0
%                              regularization, 50,51, => 1
%                              regularization, 52,53, => 2
% since the cost function is changing if using regularization
% created on Aug. 12, 2015; successful;
% modified on July,2,2016; done

function [ cost_fcn ] = Compute_Cost(para, f, g, r, h, scalar)
    
    [T,N] = size(f);
    t = linspace(0,1,T)';    
    cost_fcn = 0;
    GA = para.GA_choice;
    
    if nargin > 6
        warning('Too many intputs in Compute_Cost\n');
    end
    
    for i = 1:N
        %%% all model have this %%%
        cost_fcn = cost_fcn + ...
            SquareL2norm( t, f(:,i) - scalar(1,i)*MyGroupAction(t, g(:,1), r(:,i),GA) - h(:,1) );
    end
    cost_fcn = cost_fcn/N;
end

