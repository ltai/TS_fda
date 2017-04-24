% this is a routine for generating g and h
% created at Sep. 9, 2015 

clear; close all;

T = 201;

t = linspace(0,1,T)';

g = cos(10*pi*t);
h = 5*exp( -20*(t-0.5).^2 );          % 5*exp( -20*(t-0.5).^2 )

[approx_g,~] = L2_Shift_Legendre_Approximation(g,5,20);
[approx_h,~] = L2_Shift_Legendre_Approximation(h,1,4);