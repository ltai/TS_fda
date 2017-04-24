% this routine computes L1 norm of a function between [0,1]

function [ L1norm ] = L1norm( t,f )
    a = 0; b = 1;
    L1norm = Simpsons_Rule( abs(f),a,b );
end

