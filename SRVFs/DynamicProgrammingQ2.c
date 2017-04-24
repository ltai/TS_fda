/* input two functions, find the optimal gamma (warping function). 
 * usage: gam = DynamicProgrammincQ(q1,q2,0,0); 
 * Align q1 to q2; where q is the SRVF space 
 * 
 * commment:  this is from Wen's version, the only difference is when you get approximated gramma
 *            from the DP, then normalize from [1/n,1] to [0,1]. The Srivastava's version is
 *            doing outside DP.
 *           
 *            This problem is solving min|q_1 - q_2(r)sqrt(dr)|_L2 
 *
 * 
 * Modified: Oct.1, 2013,  add a lot of comments      
 * Check: Jan. 12, 2014,   change some variable names    
 * Check: Jan. 16, 2014    add c environment
 * Modified: June, 7, 2014 begin to change data structure, 
 *                         free to choose backward dependent set   
 *
 */

#ifndef DynamicProgrammingQ2_C
#define DynamicProgrammingQ2_C

/*===== define environment MEX_ENV or C_ENV ===== 
 * currently only finished MEX_ENV, 
 * C_ENV will be future goal, for the purpose of speed
 * */
 
#define MEX_ENV

#ifdef MEX_ENV
#include "mex.h"       /* mex.h is used for matlab */
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NNBRS 23 /* change number of grid points here!
                  choices are: 3, 7, 11, 23, 35, 63, 91*/

#if NNBRS == 3
 int Nbrs[NNBRS][2] = { 
	{ 1, 1 },{ 1, 2 },{ 2, 1 }
}

#elif NNBRS == 7
 int Nbrs[NNBRS][2] = { 
	{ 1, 1 },{ 1, 2 },{ 2, 1 },{ 2, 3 },{ 3, 2 },{ 1, 3 },
	{ 3, 1 }
};

#elif NNBRS == 11
 int Nbrs[NNBRS][2] = { 
	{ 1, 1 },{ 1, 2 },{ 2, 1 },{ 2, 3 },{ 3, 2 },{ 1, 3 },
	{ 3, 1 },{ 1, 4 },{ 3, 4 },{ 4, 3 },{ 4, 1 }
};

#elif NNBRS == 23
 int Nbrs[NNBRS][2] = { 
	{ 1, 1 },{ 1, 2 },{ 2, 1 },{ 2, 3 },{ 3, 2 },{ 1, 3 },
	{ 3, 1 },{ 1, 4 },{ 3, 4 },{ 4, 3 },{ 4, 1 },{ 1, 5 },
	{ 2, 5 },{ 3, 5 },{ 4, 5 },{ 5, 4 },{ 5, 3 },{ 5, 2 },
	{ 5, 1 },{ 1, 6 },{ 5, 6 },{ 6, 5 },{ 6, 1 }
};

#elif NNBRS == 35
  int Nbrs[NNBRS][2] = {
    { 1,  1 }, { 1,  2 }, { 1,  3 }, { 1,  4 }, { 1,  5 }, { 1,  6 }, { 1,  7 },
    { 2,  1 }, { 2,  3 }, { 2,  5 }, { 2,  7 }, 
    { 3,  1 }, { 3,  2 }, { 3,  4 }, { 3,  5 }, { 3,  7 },
    { 4,  1 }, { 4,  3 }, { 4,  5 }, { 4,  7 }, 
    { 5,  1 }, { 5,  2 }, { 5,  3 }, { 5,  4 }, { 5,  6 }, { 5,  7 }, 
    { 6,  1 }, { 6,  5 }, { 6,  7 }, 
    { 7,  1 }, { 7,  2 }, { 7,  3 }, { 7,  4 }, { 7,  5 }, { 7,  6 } 
};

#elif NNBRS == 63
int Nbrs[NNBRS][2] = {
    { 1,  1 }, { 1,  2 }, { 1,  3 }, { 1,  4 }, { 1,  5 }, { 1,  6 }, { 1,  7 }, { 1,  8 }, { 1,  9 }, { 1, 10 }, 
    { 2,  1 }, { 2,  3 }, { 2,  5 }, { 2,  7 }, { 2,  9 }, 
    { 3,  1 }, { 3,  2 }, { 3,  4 }, { 3,  5 }, { 3,  7 }, { 3,  8 }, { 3, 10 }, 
    { 4,  1 }, { 4,  3 }, { 4,  5 }, { 4,  7 }, { 4,  9 }, 
    { 5,  1 }, { 5,  2 }, { 5,  3 }, { 5,  4 }, { 5,  6 }, { 5,  7 }, { 5,  8 }, { 5,  9 }, 
    { 6,  1 }, { 6,  5 }, { 6,  7 }, 
    { 7,  1 }, { 7,  2 }, { 7,  3 }, { 7,  4 }, { 7,  5 }, { 7,  6 }, { 7,  8 }, { 7,  9 }, { 7, 10 }, 
    { 8,  1 }, { 8,  3 }, { 8,  5 }, { 8,  7 }, { 8,  9 }, 
    { 9,  1 }, { 9,  2 }, { 9,  4 }, { 9,  5 }, { 9,  7 }, { 9,  8 }, { 9, 10 }, 
    {10,  1 }, { 10,  3 }, { 10,  7 }, { 10,  9 } 
};


#elif NNBRS	== 91 
 int Nbrs[NNBRS][2] = { 
    { 1,  1}, { 1,  2}, { 1,  3}, { 1,  4}, { 1,  5}, { 1,  6},
    { 1,  7}, { 1,  8}, { 1,  9}, { 1, 10}, { 1, 11}, { 1, 12},
    { 2,  1}, { 2,  3}, { 2,  5}, { 2,  7}, { 2,  9}, { 2, 11},
    { 3,  1}, { 3,  2}, { 3,  4}, { 3,  5}, { 3,  7}, { 3,  8},
    { 3, 10}, { 3, 11}, 
    { 4,  1}, { 4,  3}, { 4,  5}, { 4,  7}, { 4,  9}, { 4, 11},
    { 5,  1}, { 5,  2}, { 5,  3}, { 5,  4}, { 5,  6}, { 5,  7},
    { 5,  8}, { 5,  9}, { 5, 11}, { 5, 12},
    { 6,  1}, { 6,  5}, { 6,  7}, { 6, 11},
    { 7,  1}, { 7,  2}, { 7,  3}, { 7,  4}, { 7,  5}, { 7,  6},
    { 7,  8}, { 7,  9}, { 7, 10}, { 7, 11}, { 7, 12},
    { 8,  1}, { 8,  3}, { 8,  5}, { 8,  7}, { 8,  9}, { 8, 11},
    { 9,  1}, { 9,  2}, { 9,  4}, { 9,  5}, { 9,  7},
    { 9,  8}, { 9, 10}, { 9, 11},
    {10,  1}, {10,  3}, {10,  7}, {10,  9}, {10, 11},
    {11,  1}, {11,  2}, {11,  3}, {11,  4}, {11,  5}, {11,  6},
    {11,  7}, {11,  8}, {11,  9}, {11, 10}, {11, 12},
    {12,  1}, {12,  5}, {12,  7}, {12, 11}
};
#endif /* NNBRS */

double CostFnDTW(const double *q1, const double *q2);
double CostFn2(const double *q1, const double *q2, const double *q2L, int k,\
                    int l, int i, int j, int n, int N, int M, double lam);
                    
void thomas(double *x, const double *a, const double *b, double *c, int n);
void spline(double *D, const double *y, int n);
void lookupspline(double *t, int *k, double dist, double len, int n);
double evalspline(double t, const double D[2], const double y[2]);

#ifdef MEX_ENV
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
/* In Matlab format:
 *    [gam,Emin] = DynamicProgrammingQ(const double *p_q1, const double *p_q2, 0, 0)
 * Input Argument:
 *     q1    double, Nx1, SRVF of f1 
 *     q2    double, Nx1, SRVF of f2
 *     lambda      scalar, coef. of smoothness panelty
 *     disp_opt    scalar, display or not
 * Output Argument:
 * 	   gamma    double, Nx1, warping path? or anything else
 *     Emin     minimized accumulated distance
 * Intermediate Argument:
 *     n     # of rows, # of functions
 *     N     # of cols, # of time grids
 *     M     = n*5*N,  
 * 	   q2L   vector , interpolate value of q2 ?
 *     D     vector, 
 *     temp  vector, 	
 *     E     is energy
 * 	   x     real wapring path in x
 *     y     real warping path in y
 *     N     size of q1, same size of q2          
 *     path  of size 2*N*N, save warping path index
 *     Eidx  save the index that has minimum value
 *     cnt   the number of real warping path
 *     (l,k)  is local path start point
 *     (i,j)  is local path end point
 * 
 * Comment:
 *     (1)variable such as "path", "E",... are indeed 2D array, but why
 *        use 1D dynamical array. Why not just use double pointer 
 *     (2) indeed, the DTW here is a simpler version, we assume the size of
 *         two input data are the same with uniformly spacing. Moreover,
 *         no adjustment window. 
 */
	int i, j, k, l, n, N, M, Eidx, Fidx, Ftmp, Fmin, Num, NumChoice, *Path, *x, *y, cnt;
	const double *q1, *q2;
	double *q2L, *D, *tmp, *E, Etmp, Emin, t, a, b, lam = 0,*gamma;
	
    if (nrhs != 4)
		mexErrMsgTxt("usage: [gam] = DynamicProgrammingQ(q1,q2,(defunct)lam,(defunct)Disp)");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("Expected double precision arguments.");

	n = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	
	if (n != mxGetM(prhs[1]) || N != mxGetN(prhs[1]))
		mexErrMsgTxt("Dimension mismatch between first and second argument.");

	if (nlhs > 1)
		mexErrMsgTxt("Expected one return.");
    
	/* temparaly not add any debug check routine */
    q1 = mxGetPr(prhs[0]);	
	q2 = mxGetPr(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);	
	/*plhs[2] = mxCreateNumericArray(1, N, mxINT);
	plhs[3] = mxCreateNumericArray(1, N, mxINT);*/
	gamma = mxGetPr(plhs[0]);

    /* initialization of spline stuff */
	M = 5*N;                                          /* change 5 to 10 ? higher accuracy */
	q2L = (double *)malloc(n*M*sizeof(double));
	D = (double *)malloc(2*N*sizeof(double));
	tmp = D + N;
    
	a = 1.0/N; /* left grid point */
	b = 1.0;   /* right grid point */
    
	for (i = 0; i < n; ++i) {
		for (j = 0; j < N; ++j) {
			tmp[j] = q2[n*j + i];
		}
		spline(D, tmp, N);

		for (j = 0; j < M; ++j) {
			/* Extrapolation 1/M < 1/N, I don't know why here */
			lookupspline(&t, &k, (j+1.0)/M - a, b - a, N);
			q2L[n*j + i] = evalspline(t, D+k, tmp+k);
		}
	}
	free(D);
	
    /* Initialization Energy and Path */
	E = (double *)calloc(N*N, sizeof(double)); /* calloc is with I.V. zero */
	Path = (int *)malloc(2*N*N*sizeof(int));

	for (i = 0; i < N; ++i) {
		E[N*i + 0] = 100;
		E[N*0 + i] = 100;
		
		Path[N*(N*0 + i) + 0] = -1;
		Path[N*(N*0 + 0) + i] = -1;
		Path[N*(N*1 + i) + 0] = -1;
		Path[N*(N*1 + 0) + i] = -1;
	}
	E[N*0 + 0] = 0;
	
    /* do the minimization at each grid point, construct the local
       distance matrix.  */
	for (j = 1; j < N; ++j) {
		for (i = 1; i < N; ++i) {
			Emin = 100000;
			Eidx = 0;
			for (Num = 0; Num < NNBRS; ++Num){
				k = i - Nbrs[Num][0];
			    l = j - Nbrs[Num][1];
				
				/* add up cost of local path */
				if (k >= 0 && l >= 0) {
					Etmp = E[N*l + k] + CostFn2(q1,q2,q2L,k,l,i,j,n,N,M,lam);
					if (Num == 0 || Etmp < Emin) {
						Emin = Etmp;
						Eidx = Num;
					}
				}
			}
			E[N*j + i] = Emin;
			Path[N*(N*0 + j) + i] = i - Nbrs[Eidx][0];
			Path[N*(N*1 + j) + i] = j - Nbrs[Eidx][1];
		}
	}
	free(E);
	free(q2L);
	(*mxGetPr(plhs[1])) = Emin; /* output mimimized dist */
	
    /* Back Tracking to find the optimal warping function */
	x = (int *)malloc(2*N*sizeof(int)); 
	y = x + N;
	x[0] = N-1;
	y[0] = N-1;

    /* find the right size of warping path */
	cnt = 1;
	while (x[cnt-1] > 0) {
		y[cnt] = Path[N*(N*0 + x[cnt-1]) + y[cnt-1]];
		x[cnt] = Path[N*(N*1 + x[cnt-1]) + y[cnt-1]];
		++cnt;
	}
	free(Path); 
	
    /* switch order, from small to big*/
	for (i = 0, j = cnt-1; i < j; ++i, --j){
		k = x[i];
		x[i] = x[j];
		x[j] = k;

		k = y[i];
		y[i] = y[j];
		y[j] = k;
	}
    
    /*x = meGetPr(plhs[2]);
    y = meGetPr(plhs[3]);*/
    
    /* find the min at each step, save index to gamma */
	for (i = 0; i < N; ++i) {
		Fmin = 100000;
		Fidx = 0;
		for (j = 0; j < cnt; ++j) {
			Ftmp = (int)fabs(i - x[j]);
			if (j == 0 || Ftmp < Fmin) {
				Fmin = Ftmp;
				Fidx = j;
			}
		}

		if (x[Fidx] == i) {
			gamma[i] = (y[Fidx]+1);
		}
		else if(x[Fidx] > i){			
			a = x[Fidx] - i;
			b = i - x[Fidx-1];
			gamma[i] = (a*(y[Fidx-1]+1) + b*(y[Fidx]+1))/(a+b);
            /* what is it doing here? */
		}
		else{
			a = i - x[Fidx];
			b = x[Fidx+1] - i;
			gamma[i] = (a*(y[Fidx+1]+1) + b*(y[Fidx]+1))/(a+b);		
		} 
		gamma[i] /= N;
	}

	/* normalize gamma again, make sure in interval [0,1], not in original version 
     make sure this is correct. I expect gamma[0]=0 in the last line, not here
     what's wrong? */
    
    a = gamma[0];
    b = gamma[N - 1] - gamma[0];
    for(i = 0; i < N; i++){
        gamma[i] = (gamma[i] - a) / b;
    }
	free(x);
}
#endif

#ifdef C_ENV{
/* NOT finished!!!!!!!!!!!!!!! Jan. 16, 2014
 * The material is nothing different from the MEX_ENV version
 * */
 
   	int i, j, k, l, n, N, M, Eidx, Fidx, Ftmp, Fmin, Num, NumChoice, *Path, *x, *y, cnt;
	const double *q1, *q2;
	double *q2L, *D, *tmp, *E, Etmp, Emin, t, a, b, lam = 0,*gamma;
	
	/* temparaly not add any debug check routine */
	n = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
    q1 = mxGetPr(prhs[0]);	
	q2 = mxGetPr(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);	
	gamma = mxGetPr(plhs[0]);
    
    /* initialization of spline stuff */
	M = 5*N;
	q2L = (double *)malloc(n*M*sizeof(double));
	D = (double *)malloc(2*N*sizeof(double));
	tmp = D + N;
	a = 1.0/N;
	b = 1.0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < N; ++j) {
			tmp[j] = q2[n*j + i];
		}
		spline(D, tmp, N);

		for (j = 0; j < M; ++j) {
			/* XXX: Extrapolation 1/M < 1/N */
			lookupspline(&t, &k, (j+1.0)/M - a, b - a, N);
			q2L[n*j + i] = evalspline(t, D+k, tmp+k);
		}
	}
	free(D);
	
    /* Initialization Energy and Path */
	E = (double *)calloc(N*N, sizeof(double)); 
	/* same as malloc, but with initial value zero */
	Path = (int *)malloc(2*N*N*sizeof(int));

	for (i = 0; i < N; ++i) {
		E[N*i + 0] = 1;
		E[N*0 + i] = 1;
		
		Path[N*(N*0 + i) + 0] = -1;
		Path[N*(N*0 + 0) + i] = -1;
		
		Path[N*(N*1 + i) + 0] = -1;
		Path[N*(N*1 + 0) + i] = -1;
	}
	E[N*0 + 0] = 0;
	
    /* do the minimization at each grid point, construct the local
     * distance matrix. This part should expect DTW 
     * and the Elastic Metric are the same 
     */
	for (j = 1; j < N; ++j) {
		for (i = 1; i < N; ++i) {
			Emin = 100000;
			Eidx = 0;
			for (Num = 0; Num < NNBRS; ++Num){
				k = i - Nbrs[Num][0];
			    l = j - Nbrs[Num][1];
			    
				/* choose right dependent set */
				/*if(choice==0){
					k = i - NDTW[Num][0];
					l = j - NDTW[Num][1];					
				}
				else if(choice==1){
					k = i - Nbrs[Num][0];
					l = j - Nbrs[Num][1];
				}
				else printf("Wrong in DP code");*/
				
				/* add up cost of local path */
				if (k >= 0 && l >= 0) {
					Etmp = E[N*l + k] + CostFn2(q1,q2,q2L,k,l,i,j,n,N,M,lam);
					if (Num == 0 || Etmp < Emin) {
						Emin = Etmp;
						Eidx = Num;
					}
				}
			}
			E[N*j + i] = Emin;
			Path[N*(N*0 + j) + i] = i - Nbrs[Eidx][0];
			Path[N*(N*1 + j) + i] = j - Nbrs[Eidx][1];
		}
	}
	free(E);
	free(q2L);
	
	/* output mimimized dist as r.h.s */
	(*mxGetPr(plhs[1])) = Emin;
	
    /* Back Tracking to find the optimal warping function */
	x = (int *)malloc(2*N*sizeof(int));
	y = x + N;
	x[0] = N-1;
	y[0] = N-1;

	cnt = 1;
	while (x[cnt-1] > 0) {
		y[cnt] = Path[N*(N*0 + x[cnt-1]) + y[cnt-1]];
		x[cnt] = Path[N*(N*1 + x[cnt-1]) + y[cnt-1]];
		++cnt;
	}
	free(Path);
	
	/* right size of warping path? */
	for (i = 0, j = cnt-1; i < j; ++i, --j){
		k = x[i];
		x[i] = x[j];
		x[j] = k;

		k = y[i];
		y[i] = y[j];
		y[j] = k;
	}
    
    /* find the min at each step, save index to gamma */
	for (i = 0; i < N; ++i) {
		Fmin = 100000;
		Fidx = 0;
		for (j = 0; j < cnt; ++j) {
			Ftmp = (int)fabs(i - x[j]);
			if (j == 0 || Ftmp < Fmin) {
				Fmin = Ftmp;
				Fidx = j;
			}
		}

		if (x[Fidx] == i) {
			gamma[i] = (y[Fidx]+1);
		}
		else if(x[Fidx] > i){			
			a = x[Fidx] - i;
			b = i - x[Fidx-1];
			gamma[i] = (a*(y[Fidx-1]+1) + b*(y[Fidx]+1))/(a+b);
		}
		else{
			a = i - x[Fidx];
			b = x[Fidx+1] - i;
			gamma[i] = (a*(y[Fidx+1]+1) + b*(y[Fidx]+1))/(a+b);		
		}
		gamma[i] /= N;
	}

	/* normalize gamma again, make sure in interval [0,1] */
    a = gamma[0];
    b = gamma[N - 1] - gamma[0];
    for(i = 0; i < N; i++){
        gamma[i] = (gamma[i] - a) / b;
    }
	free(x);    
} 
#endif

double CostFn2(const double *q1, const double *q2, const double *q2L, int k,\
                    int l, int i, int j, int n, int N, int M, double lam){
/* The cost function is the most important part of DP
 * It already incorporate the SRVF concept in the coding. Watch out
 * 
 * Input
 *     q1,q1 are input magnitude value
 *     q2L  ??
 *     i,j,k,l   index for time
 * Output
 *     CostFn2  cost of path from (k,l) to (i,j)
 * Intermediate 
 *     m        slope from (k,l) to (i,j), which is derivative of gamma
 *     sqrtm    square root of derivative of gamma 
 *     ip       the integral part of a double from from modf
 *     fp       ??
 *     idx      warping index?
 *     x        ??
 *     d        ??
 *     y        line through (k,l) and (i,j)
 * */
	double m = (j-l)/(double)(i-k);
	double sqrtm = sqrt(m), E = 0, y, tmp, ip, fp;
	int x, idx, d;

	/* by assumption, gamma is linear in this interval */
	for (x = k; x <= i; ++x)
	{
		/* find the index for computing q2(r) */
		y = (x-k)*m + l + 1;
		fp = modf(y*M/N, &ip); /*" modf" is to seperate integer and fraction */
		idx = (int)(ip + (fp >= 0.5)) - 1;
        
		for (d = 0; d < n; ++d) 
		{
			tmp = q1[n*x + d] - sqrtm*q2L[n*idx + d];
			E += tmp*tmp;
		}
	}
	return E/N;
}

void thomas(double *x, const double *a, const double *b, double *c, int n){
/* Thomas Algorithm to solve a tri-diagonal matrix     
 * Input
 *     a,b  one is diagnal, the other is sup(sub) of tri-diag.  
 *     c    r.h.s. of the tri-diagonal matrix
 *     n    is the size of tri-diagonal matrix
 * Output
 *     x    the solution of tri-diagonal matrix
 * */	     
	double tmp;
	int i;

	c[0] /= b[0];
	x[0] /= b[0];

	for (i = 1; i < n; ++i) {
		tmp = 1/(b[i] - c[i-1] * a[i]);
		c[i] *= tmp;
		x[i] = (x[i] - x[i-1] * a[i])*tmp;
	}

	for (i = n-2; i >= 0; --i) {
		x[i] -= c[i]*x[i+1];
	}
}

void spline(double *D, const double *y, int n){
/* Input
 * 	  y, function value
 *    n, size of y 
 * Output
 *    D, coefficient induced by spline, solving from a tri-diagonal matrix
 * Intermediate
 *    a,b,c  the tri-diagonal matrix is s.p.d. since equally space grid,
 *           one is diag, the other is sup and sub, the rest is r.h.s. 
 *           so far I can't tell which is which
 * */
	int i;
	double *a, *b, *c;

	a = (double *)malloc(3*n*sizeof(double));
	b = a + n;
	c = b + n;

	if (n < 4) {
		a[0] = 0; 
		b[0] = 2;
		c[0] = 1;
		D[0] = 3*(y[1]-y[0]);

		a[n-1] = 1;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 3*(y[n-1]-y[n-2]);
	}
	else {
		a[0] = 0;
		b[0] = 2;
		c[0] = 4;
		D[0] = -5*y[0] + 4*y[1] + y[2];

		a[n-1] = 4;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 5*y[n-1] - 4*y[n-2] - y[n-3];
	}

	for (i = 1; i < n-1; ++i) {
		a[i] = 1;
		b[i] = 4;
		c[i] = 1;
		D[i] = 3*(y[i+1]-y[i-1]);
	}

	thomas(D, a, b, c, n);
	free(a);
}

void lookupspline(double *t, int *k, double dist, double len, int n){
/* Input
 * 	  dist
 *    len
 *    n
 * Output
 *    t
 *    k
 * */
	*t = (n-1)*dist/len;
	*k = (int)floor(*t);

	*k = (*k > 0)*(*k);
	*k += (*k > n-2)*(n-2-*k);

	*t -= *k;
}

double evalspline(double t, const double D[2], const double y[2]){
/* Input
 *     D,y   coefficients
 *     t     time you want to evaluate
 * Output
 *     evalspline
 *  */	
	double c[4];

	c[0] = y[0];
	c[1] = D[0];
	c[2] = 3*(y[1]-y[0])-2*D[0]-D[1];
	c[3] = 2*(y[0]-y[1])+D[0]+D[1];

	return t*(t*(t*c[3] + c[2]) + c[1]) + c[0];
}

#endif /* DynamicProgrammingQ2_C */
