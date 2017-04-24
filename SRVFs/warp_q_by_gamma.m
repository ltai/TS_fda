% This routine computes the warping of q (SRVF)
% another version for warp_q_by_gamma
% use the mathematical equivalence expression:
%                             (q,r) = Q(for)
% although I have test that the quivalence is indeed not hold,
% but if all the problem solving by this way. The numerical resulf if fine.
%
% created at June, 1, 2014; has tested. 
% look at Aug. 12, 2015; I think this one can be deleted since having my
%                         own routin "MyGroupAction"

function [ warp_q ] = warp_q_by_gamma( t, q, gamma )
% Input Argument:
%   t       of size Tx1
%   q       of size Tx1
%   gam     of size Tx1
% Output Argument:
%   warp_q  of size Tx1, warped f
    
    %[T,N] = size(q);
    % step 1: given q(s), compute it's corresponding f(t)
    f = cumtrapz(t, q.*abs(q));
    
    % step 2: compute for
    warp_f = warp_f_by_gamma(t,f,gamma);
       
    % step 3: compute SRVF of (for)
    warp_q = SRVF(t,warp_f); 

end
