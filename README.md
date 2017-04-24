# TVPSE
Trend and Variable-Phase Seasonality Estimation from Functional Data.
TVPSE contains all the codes in my PhD dissertation.
Programming language used were Matlab, R, and C++.

The routine "DynamicProgrammingQ.c" under sub-directory "SRVFs"
require recompile under your computer.

Step 1: open your Matlab
Step 2: Make sure current folder contains files "DynamicProgrammingQ.c".
Step 2: type mex DynamicProgrammingQ.c
Step 3: Matlab will produce "DynamicProgrammingQ.mexglx" which can be called by Matlab directly.

Further information about how to use Matlab to call c++ can be found at 
https://www.mathworks.com/help/matlab/ref/mex.html?requestedDomain=www.mathworks.com
