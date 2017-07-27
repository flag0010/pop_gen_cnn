
void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)());
void mnbrak2(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double , void *  ), double beta, void * p);
double d_trapzd(double (*func)(double), double a, double b, int n);
double d_qtrap(double (*func)(double), double a, double b);
double d_qsimp(double (*func)(double), double a, double b);
double d_qromb(double (*func)(double), double a, double b);
double midpnt2(double (*func)(double, void *), double a, double b, int n, void * p);
double d_qromb2(double (*func)(double, void *), double a, double b, void * p);
void   d_polint(double xa[], double ya[], int n, double x, 
		double *y, double *dy);
//double f(double x);
double factorial(int x);
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));
double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	     double *xmin);
double f1dim(double x);
double xchoosey(int n, int k);
double multinom(int x, int n, ...);

void nrerror(char error_text[]);
double	logfact(int num);
double	ExactProb(int nn11, int nn12, int nn21, int nn22);
double 	FsXctTst(int m11, int m12, int m21, int m22);

//void powellMin(double p[], double **xi, int n, double ftol, int *iter,
//  double *fret, double (*func)(double []))
