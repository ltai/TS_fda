% This routine compute the SRVF of input f
% create at May, 10, 2014
% modified at June, 11, 2014; gradient changes to MyGradient
% modified at June, 23, 2014, add another SRVF implementation

function [ q ] = SRVF(t, f)
%  Input   
%      t  of size Tx1
%      f  of size Tx1
%  Output 
%      q  of size Tx1, the SRVF of f
    
    binsize = mean(diff(t));
    
    [T, N] = size(f);
    %binsize = 1/(T-1);
    q = zeros(T,N);
    dev_f = zeros(T,N);
    
    for i = 1:N
        % (1) lower accuracy method
        % q(:,i) = gradient(f(:,i), binsize)./sqrt(abs(gradient(f(:,i), binsize))+eps);
        
        % (2) higher accuracy method
        % dev_f(:,i) = MyGradient(f(:,i),t);    
        % q(:,i) = dev_f(:,i)./sqrt(abs( dev_f(:,i))+eps); 
        % q(:,i) = dev_f(:,i)./sqrt(abs( dev_f(:,i))); 
        % eps for avoiding zero in the bottom
        
        % (3) another representation
        dev_f(:,i) = gradient(f(:,i),binsize);
        %dev_f(:,i) = MyGradient(f(:,i),t);
        
        q(:,i) = sign( dev_f(:,i) ).*sqrt(abs( dev_f(:,i) ));
    end
    
end

