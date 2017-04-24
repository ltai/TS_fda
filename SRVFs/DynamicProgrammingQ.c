/* This is an old vertion of DynamicProgrammingQ.c, obtained from
 * Funtional data analysis. Oct. 2012.
 * All I have done is created a lot of comments.
 * 
 * This code is solving min|q_1 -  q_2(r)sqrt(dr)|_L2
 * 
 */

#include <math.h>
#include <stdlib.h>
#include "mex.h"

/*#define NNBRS	23
const int Nbrs[NNBRS][2]={ 
	{ 1, 1 }, 
	{ 1, 2 },{ 2, 1 }, 
	{ 1, 3 },{ 2, 3 },{ 3, 2 },{ 3, 1 },
	{ 1, 4 },{ 3, 4 },{ 4, 3 },{ 4, 1 },
	{ 1, 5 },{ 2, 5 },{ 3, 5 },{ 4, 5 },
	{ 5, 4 },{ 5, 3 },{ 5, 2 },{ 5, 1 },
	{ 1, 6 },{ 5, 6 },{ 6, 5 },{ 6, 1 }
};*/

 #define NNBRS	91
const int Nbrs[NNBRS][2] = { 
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

double CostFn2(const double *q1, const double *q2, const double *q2L, \
                 int k, int l, int i, int j, int n, int N, int M, double lam);
void thomas(double *x, const double *a, const double *b, double *c, int n);
void spline(double *D, const double *y, int n);
void lookupspline(double *t, int *k, double dist, double len, int n);
double evalspline(double t, const double D[2], const double y[2]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
/*	
 * Input Argument:
 *     q1    double, Nx1, SRVF of f1
 *     q2    double, Nx1, SRVF of f2
 *     lam   constant, smoothness panalty?
 *     diap  an index to display or not, indeed, here no use
 * Output Argument:
 * 	   yy    double, Nx1, warping path? or anything else
 * 
 * Intermediate Argument:
 * 	   q2L   vector, interpolate value of q2 ?
 *     D     vector, 
 *     temp  vector, 	
 *     E     is energy?
 * 	   x     real wapring path in x
 *     y     real warping path in y
 *     N     size of q1, same size of q2      
 *     M     equal 5*N  why?????     
 *     path  of size 2*N*N, save warping path index
 *     Eidx  save the index that has minimum value
 *     cnt   the number of real warping path
 * 
 * Comment:
 *     variable such as "path", "E",... are indeed 2D array, but why
 *     auther use 1D dynamical array. Why not just use double pointer 
 */
	
	int i, j, k, l, n, M, N, Eidx, Fidx, Ftmp, Fmin, Num, *Path, *x, *y, cnt;
	const double *q1, *q2;
	double *q2L, *yy, *D, *tmp, *E, Etmp, Emin, t, a, b, lam;

	if (nrhs != 4)
		mexErrMsgTxt("usage: [gam] = DynamicProgrammingQ(q1,q2,(defunct)lam,(defunct)Disp)");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("Expected double precision arguments.");

	if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2)
		mexErrMsgTxt("First two arguments expected to be two dimensional.");

	n = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	
	if (n != mxGetM(prhs[1]) || N != mxGetN(prhs[1]))
		mexErrMsgTxt("Dimension mismatch between first and second argument.");

	if (nlhs > 1)
		mexErrMsgTxt("Expected one return.");

	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);

	yy = mxGetPr(plhs[0]);
	q1 = mxGetPr(prhs[0]);
	q2 = mxGetPr(prhs[1]);
    lam = *mxGetPr(prhs[2]);
    
	M = 50*N;  /* M = 5*N */
	q2L = malloc(n*M*sizeof(double));
	D = malloc(2*N*sizeof(double));
	tmp = D + N;

	a = 1.0/N;
	b = 1.0;

	for (i = 0; i < n; ++i){
		for (j = 0; j < N; ++j){
			tmp[j] = q2[n*j + i];
		}
		spline(D, tmp, N);
		for (j = 0; j < M; ++j) {
			/* Extrapolation 1/M < 1/N */
			lookupspline(&t, &k, (j+1.0)/M - a, b - a, N);
			q2L[n*j + i] = evalspline(t, D+k, tmp+k);
		}
	}
	free(D);
	
	/* Initialization Energy and Path */
	E = calloc(N*N, sizeof(double));
	Path = malloc(2*N*N*sizeof(int));

	for (i = 0; i < N; ++i) {
		E[N*i + 0] = 1000000;
		E[N*0 + i] = 1000000;
		Path[N*(N*0 + i) + 0] = -1;
		Path[N*(N*0 + 0) + i] = -1;
		Path[N*(N*1 + i) + 0] = -1;
		Path[N*(N*1 + 0) + i] = -1;
	}
	E[N*0 + 0] = 0;
	
	/* do the minimization at each grid point */
	for (j = 1; j < N; ++j) {
		for (i = 1; i < N; ++i) {
			Emin = 100000;
			Eidx = 0;
			for (Num = 0; Num < NNBRS; ++Num) {
				k = i - Nbrs[Num][0];
				l = j - Nbrs[Num][1];

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

	/* x, y assumed to be at most length N 
	 * convet temp path to real path ?
	 * make the size of path right
	 * */
	x = malloc(2*N*sizeof(int));
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
	
	/* Trace back to find the optimal warping function */
	for (i = 0, j = cnt-1; i < j; ++i, --j) {
		k = x[i];
		x[i] = x[j];
		x[j] = k;

		k = y[i];
		y[i] = y[j];
		y[j] = k;
	}

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
			yy[i] = (y[Fidx]+1);
		}
		else {
			if (x[Fidx] > i) {
				a = x[Fidx] - i;
				b = i - x[Fidx-1];
				yy[i] = (a*(y[Fidx-1]+1) + b*(y[Fidx]+1))/(a+b);
			}
			else {
				a = i - x[Fidx];
				b = x[Fidx+1] - i;
				yy[i] = (a*(y[Fidx+1]+1) + b*(y[Fidx]+1))/(a+b);
			}
		}
		yy[i] /= N;
	}
	free(x);
}

double CostFn2(const double *q1, const double *q2, const double *q2L,\
           int k, int l, int i, int j, int n, int N, int M, double lam){
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
 *     idx      ??
 *     x        ??
 *     d        ??
 *     y        line through (k,l) and (i,j)
 * */
			   
	double m, sqrtm, E, y, tmp, ip, fp;
	int x, idx, d;
	
	 m = (double)(j-l)/(double)(i-k);
	 sqrtm = sqrt(m);
	 E=0;

	for (x = k; x <= i; ++x){
		y = (x-k)*m + l + 1;
		fp = modf(y*M/N, &ip);     
		/*" modf" is to seperate integer and fraction */
		idx = (int)(ip + (fp >= 0.5)) - 1;
		/* is here doing approximation !!? */

		for (d = 0; d < n; ++d) {
			tmp = q1[n*x + d] - sqrtm*q2L[n*idx + d];
            E += tmp*tmp + lam*(1-sqrtm)*(1-sqrtm);
		}
	}
	return E/N;
}

/* Thomas Algorithm to solve a tri-diagonal matrix */
void thomas(double *x, const double *a, const double *b, double *c, int n){
/* Input
 *     a,c  one is diagnal, the other is sup(sub) of tri-diag.  
 *     b    r.h.s. of the tri-diagonal matrix
 *     n    is the size of tri-diagonal matrix
 * Output
 *     x    0~(n-1), the solution of tri-diagonal matrix, coef of S''
 * */	
	double tmp;
	int i;

	c[0] /= b[0];
	x[0] /= b[0];

	for (i = 1; i < n; ++i)
	{
		tmp = 1/(b[i] - c[i-1] * a[i]);
		c[i] *= tmp;
		x[i] = (x[i] - x[i-1] * a[i])*tmp;
	}
	
	for (i = n-2; i >= 0; --i) 
		x[i] -= c[i]*x[i+1];
}

void spline2(double *D, const double *y, int n){
/* Input
 * 	  y, function values
 *    n, size of y
 * Output
 *    D, coefficient induced by spline, solving a tri-diagonal matrix
 * Intermediate
 *    a,c       sub-diag and sup-diag entry
 *    b         diagonal entry
 * */
	int i;
	double *a, *b, *c;

	a = malloc(3*n*sizeof(double));
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
		/* here is the boundary condition of spline,
		 * read it!! */
		 
		a[0] = 0; 		
        b[0] = 0;
		c[0] = 1; 		
        D[0] = 0;

		a[n-1] = 1; 		
        b[n-1] = 0;
		c[n-1] = 0; 		
        D[n-1] = 0;
	}

	for (i = 1; i < n-1; ++i) {
		a[i] = 1; 		
        b[i] = 4;
		c[i] = 1; 		
        D[i] = 6*(y[i+1]-y[i-1]);
	}
	thomas(D, a, b, c, n);
	free(a);
}

/* Input
 * 	  y, function values
 *    n, size of y
 * Output
 *    D, coefficient induced by spline, solving a tri-diagonal matrix
 * Intermediate
 *    a,c       sub-diag and sup-diag entry
 *    b         diagonal entry
 * */
void spline(double *D, const double *y, int n){

	int i;
	double *a, *b, *c;

	a = malloc(3*n*sizeof(double));
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
 *    k     index of spline that you are going to evaluate
 * */
	*t = (n-1)*dist/len;
	*k = (int)floor(*t);
	
	/* what do you mean here! */
	*k = (*k > 0)*(*k);
	*k += (*k > n-2)*(n-2-*k);

	*t -= *k;
}

double evalspline(double t, const double D[2], const double y[2]){
/* Input
 *     D, S'' coef
 *     y, function value 
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
