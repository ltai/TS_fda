% this program performs normalization of observed data
% one is shift to [0,1]
%
% the other is shift to integrate(data) = 0 on [0,1]
% created on March 30, 2015

close all;
normalize = 2;

if normalize == 1      % nonlinear normalized to [0,1] 
    minf = min(f);
    maxf = max(f);

    for i=1:14
        ff(:,i) = ( f(:,i) - minf(i) )/(maxf(i)-minf(i)); 
    end
    plot(t,ff);
    
elseif normalize == 2  % linear normalized to int(f)=0
    
    [~,N] = size(f);
    
    for i = 1:N
        
        f(:,i) = f(:,i) - trapz(t,f(:,i));
        
    end
    clear N normalize i;
    plot(t,f);
end
