# TS_fda
TS_fda stands for Trend and Variable-Phase Seasonality Estimation from Functional Data.
TS_fda contains all the codes in my PhD dissertation and the paper (http://arxiv.org/abs/1704.07358). 
Programming language used were Matlab, R, and C++.

"main_function" is the code for trend and seasonality estimation given a subspace H.
Figure results are saved in a new directory "html" <br>

"Bootstrap_Experiment" serves as two purposes: <br>
(1) When set B=1, it is a subplace H selection among Fourier, sin, cos, and Legendre, with basis eleements from l1 to l2.<br>
(2) When donig Bootstrap analysis, B is B=100 or 200 depends on users' need.

------------------------------------ important ---------------------------------------<br>
The routine "DynamicProgrammingQ.c" under sub-directory "SRVFs" requires recompiling under your computer.

Step 1: open your Matlab <br>
Step 2: Make sure current folder contains the file "DynamicProgrammingQ.c" <br>
Step 3: type mex DynamicProgrammingQ.c <br>
Step 4: Matlab will produce "DynamicProgrammingQ.mexglx" which can be used by Matlab directly. <br>

Further information about how to use Matlab to call c++ can be found at 
https://www.mathworks.com/help/matlab/ref/mex.html?requestedDomain=www.mathworks.com
