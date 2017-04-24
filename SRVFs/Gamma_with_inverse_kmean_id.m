%  This routine gives a set of gamma whose inverse gamma has
%  Karcher mean r_id
%  created at April, 29, 2014
%  continued on Sep. 1, 2015; my version works!! here after, use idx = 2.

function [ imposed_gam, mu, muI ] = Gamma_with_inverse_kmean_id( gam )
%  input    r     of size TxN; M is time; N is # of gammas
%  output   gam   of size TxN;

    [T,N] = size(gam);
    t = linspace(0,1,T)';
    
    idx = 2; % 1 is original paper version; 2 is my own version 
    
    if idx == 1
        imposed_gam = zeros(T,N);
        gamI = zeros(T,N);

        for i = 1:N
            gamI(:,i) = invertGamma(gam(:,i));
        end

        mu = Karcher_Mean_Gamma(gamI);
        muI = invertGamma(mu);

        for i = 1:N
            gamI(:,i) = interp1(t,gamI(:,i),muI); 
            if( isnan(gamI(1,i))|| isnan(gamI(end,i)))
                disp('wrong in subroutine Gamma_with_inverse_keman_id');
            end
        end

        for i = 1:N
            imposed_gam(:,i) = invertGamma( gamI(:,i) );
        end
        
    elseif idx == 2
        imposed_gam = zeros(T,N);
        gamI = zeros(T,N);
        
        for i = 1:N
            gamI(:,i) = invertGamma(gam(:,i));
        end
        
        mu = Karcher_Mean_Gamma(gamI);
        muI = invertGamma(mu);
        
        for i = 1:N
           imposed_gam(:,i) = interp1(t,mu,gam(:,i)); 
           if isnan(imposed_gam(end,i))
               imposed_gam(end,i) = 1;
           end 
        end
        
        for i = 1:N
            if( isnan(imposed_gam(1,i))|| isnan(imposed_gam(end,i)))
                disp('NaN in Gamma_with_inverse_keman_id');
            end
        end
        
    end
end

