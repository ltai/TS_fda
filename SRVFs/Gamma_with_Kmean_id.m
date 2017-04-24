%  This routine gives a set of gamma whose Karcher mean is r_id
%  refer to Srivastava, 2011, Lemma 4; not going to phi = r^2 space
%  An approximation, or alternative of Algorithm 1
%  created at June, 2, 2014; to be rewrited
%

function [ gam ] = Gamma_with_Kmean_id( r )
%  Input    r     of size TxN; M is time; N is # of gammas
%  Output   gam   of size TxN;

    [T,N] = size(r);
    t = linspace(0,1,T)';
    gam = zeros(T,N);
    
    mean_gam = mean(r')';
    mean_gamI = invertGamma(mean_gam);
    
    for i=1:N
        gam(:,i) = interp1(t, r(:,i), mean_gamI);
        gam(:,i) = gam(:,i)/gam(T,i);
    end
    
end


