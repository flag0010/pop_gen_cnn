/* numerical stuff */

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "numerical.h"
#include <gsl/gsl_sf.h>
//#include "nrutil.h"

/* Numerical Recipes Stuff */

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY  1.0e-20
#define MAX(a,b) ((a) > (b)  ? (a) : (b))
#define MIN(a,b) ((a) < (b)  ? (a) : (b))
#define SIGN(a,b)  ((b)  > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define FUNC(x) ((*func)(x))
#define FUNC2(x,p) ((*func)(x,p))
#define MIN_P_VALUE 1.0e-200
#define EPS 1.0e-5
#define JMAX 50


/* this is from numerical recipes; brackets a minimum with ax, bx, cx */
void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)()){
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if(*fb > *fa){
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = *bx + GOLD * (*bx-*ax);
	*fc=(*func)(*cx);
	while(*fb > *fc){
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u = *bx-((*bx-*cx)*q-(*bx-*ax)*r) /
		    (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = *bx + GLIMIT * (*cx-*bx);
		if((*bx-u) * (u-*cx) > 0.0){
			fu=(*func)(u);
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
   	 			}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return;
				}
			u = *cx + GOLD * (*cx-*bx);
			fu = (*func)(u);
			}
		else if((*cx-u)*(u-ulim) > 0.0){
			fu=(*func)(u);
			if(fu < *fc){
				SHFT(*bx,*cx,u,*cx + GOLD * (*cx - *bx));
				SHFT(*fb,*fc,fu,(*func)(u));
			}
		}
		else if((u-ulim)*(ulim-*cx) >= 0.0){
			u=ulim;
			fu=(*func)(u);
		}
		else{
			u = *cx+GOLD*(*cx-*bx);
			fu=(*func)(u);
			}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
		}
	return;
}

/* this is from numerical recipes; brackets a minimum with ax, bx, cx */
void mnbrak2(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double, void * ), double beta, void * p){
	double ulim,u,r,q,fu,dum;
	

	*fa=(*func)(*ax, p);
	*fb=(*func)(*bx, p);
	if(*fb > *fa){
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = *bx + GOLD * (*bx-*ax);
	*fc=(*func)(*cx, p);
	while(*fb > *fc){
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u = *bx-((*bx-*cx)*q-(*bx-*ax)*r) /
		    (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = *bx + GLIMIT * (*cx-*bx);
		if((*bx-u) * (u-*cx) > 0.0){
			fu=(*func)(u, p);
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
   	 			}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return;
				}
			u = *cx + GOLD * (*cx-*bx);
			fu = (*func)(u, p);
			}
		else if((*cx-u)*(u-ulim) > 0.0){
			fu=(*func)(u, p);
			if(fu < *fc){
				SHFT(*bx,*cx,u,*cx + GOLD * (*cx - *bx));
				SHFT(*fb,*fc,fu,(*func)(u, p));
			}
		}
		else if((u-ulim)*(ulim-*cx) >= 0.0){
			u=ulim;
			fu=(*func)(u, p);
		}
		else{
			u = *cx+GOLD*(*cx-*bx);
			fu=(*func)(u, p);
			}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
		}
	return;
}


/* Numerical Integration Routines */


#define FUNC(x) ((*func)(x))
double d_trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

/* specific version for GSL function */
double d_trapzd2(double (*func)(double, void *), double a, double b, int n, void * p)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
	  return (s=0.5*(b-a)*(FUNC2(a,p)+FUNC2(b,p)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC2(x,p);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

double d_qtrap(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double s,olds;
	int j;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=d_trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}

/* specific version for GSL function */
double d_qtrap2(double (*func)(double, void *), double a, double b, void * p){
  double s,olds;
  int j;

  olds = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    s=d_trapzd2(func,a,b,j,p);
    if (fabs(s-olds) < EPS*fabs(olds)) return s;
    olds=s;
  }
  nrerror("Too many steps in routine qtrap");
  return 0.0;
}

double d_qsimp(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	int j;
	double s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=d_trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}

#define JMAXP (JMAX+1)
#define K 5
double d_qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double ss,dss;
	double s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=d_trapzd(func,a,b,j);
		if (j >= K) {
			d_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}

double d_qromb2(double (*func)(double, void *), double a, double b, void * p){

  double ss,dss;
  double s[JMAXP+1],h[JMAXP+1];
  int j;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=midpnt2(func,a,b,j,p);
    if (j >= K) {
      d_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef JMAXP
#undef K

#define NRANSI
#include "nrutil.h"
void d_polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}
#undef NRANSI



/* from chp. 4 of numerical recipes */
double midpnt(double (*func)(double), double a, double b, int n){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC(0.5*(a+b)));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}

/* from chp. 4 of numerical recipes, specific for GSL functions */
double midpnt2(double (*func)(double, void *), double a, double b, int n, void * p){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC2(0.5*(a+b),p));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC2(x,p);
      x += ddel;
      sum += FUNC2(x,p);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}


/* factorials */

double xchoosey(int n, int k){
  if(n<k) return(0);
  if(n == k || k == 0) return(1);
  return(gsl_sf_choose(n,k));
}

double multinom(int x, int n, ...) {
  int i, y;
  va_list ap;
  double ret_val=0.0;
  
  va_start(ap, n);
  for (i=0; i<n; i++) {
    y = va_arg(ap, int);
    if (x < y) return 0.0;
    ret_val *= xchoosey(x, y);
    x -= y;
    if (x < 0) return 0.0;
  }
  va_end(ap);
  return exp(ret_val);
}


//from Chuck Langley
/* +++++++++++  chl   19JAN2010  ++++++++++++++++++++++++++++++++++++++++++  */

double FsXctTst(int m11, int m12, int m21, int m22) {
	double      onetail, twotail;
	double		ep[(m11+m12+m21+m22 + 1)];
	double		cov[(m11+m12+m21+m22 + 1)];	
	int			nt11,nt12,nt21,nt22 ;
	int			i ;
	if( ((m11+m21)<1) || ((m12+m22)<1) )   return(-99.0);
	nt11 = 0 ;
	do
	{
		nt12 = m11 + m12 - nt11 ;
		nt21 = m11 + m21 - nt11 ;
		nt22 = m11 + m12 + m21 + m22 - nt11 - nt12 - nt21 ;
		
		ep[nt11] = 0.0 ;
		cov[nt11] = 0.0 ;
		
		if ((nt12 < 0) || (nt21 < 0) || (nt22 < 0))
			nt11++ ;
		else
		{
			ep[nt11] = ExactProb(nt11,nt12,nt21,nt22) ;
			cov[nt11] = nt11*nt22 - nt12*nt21 ;
			
			/*	fprintf(stderr,"\nep[%2d] = %e\tcov[%2d] = %e", nt11,ep[nt11],nt11,cov[nt11]) ;	*/
			nt11++ ;
		}
	}	while (nt11 <= (MIN( (m11 + m12),(m11 + m21))) ) ;

	if(cov[m11] <= 0.0){
		onetail = 0.0 ;
		
		for(i=0;i<=m11;i++)	onetail += ep[i] ;
		
		twotail = onetail ;
		i = nt11 - 1 ;
		if(ep[i] <= (ep[m11] + 0.0001))
		{
			do
			{
				twotail += ep[i] ;
				i-- ;
			}	while ((ep[i] <= (ep[m11] + 0.0001)) && (i > m11)) ;
		}
	}
	else
	{
		onetail = 0.0 ;
		
		for(i=m11;i<nt11;i++)	onetail += ep[i] ;
		
		twotail = onetail ;
		i = 0 ;
		if(ep[0] <= (ep[m11] + 0.0001))
		{
			do
			{
				twotail += ep[i] ;
				i++ ;
			}	while ((ep[i] <= (ep[m11] + 0.0001)) && (i < m11)) ;
		}
	}
	if(cov[m11] == 0.0)	onetail = 99 ;
	
	/*	fprintf(stderr,"\n onetail = %e \t twotail = %e ", *onetail, *twotail) ;		*/
	if(twotail < MIN_P_VALUE)   printf("\nFsXctTst(%d,%d,%d,%d) = %lf", m11, m12, m21, m22, twotail);
	return(twotail) ;
}


double	ExactProb(int nn11, int nn12, int nn21, int nn22){
	double	xp = 0.0;
	
	xp  = logfact(nn11+nn12) ;
	xp += logfact(nn11+nn21) ;
	xp += logfact(nn12+nn22) ;
	xp += logfact(nn21+nn22) ;
	xp -= logfact(nn11+nn12+nn21+nn22) ;
	xp -= logfact(nn11) ;
	xp -= logfact(nn12) ;
	xp -= logfact(nn21) ;
	xp -= logfact(nn22) ;
	
	return(exp(xp)) ;
}

double	logfact(int num){
	double	sumlog = 0.0 ;
	
	if((num == 0) || (num == 1))	return(0.0) ;
	else
	{
		do
		{
			sumlog += log(((double)num)) ;
			num-- ;
		} while (num > 1) ;
		
		return(sumlog) ;
	}
    return(-99.99) ;
}


/* #define ITMAX 200 */
/* #define NRANSI */

/* void powellMin(double p[], double **xi, int n, double ftol, int *iter, double *fret, */
/* 	double (*func)(double [])) */
/* { */
/* 	void linmin(double p[], double xi[], int n, double *fret, */
/* 		float (*func)(double [])); */
/* 	int i,ibig,j; */
/* 	float del,fp,fptt,t,*pt,*ptt,*xit; */

/* 	pt=vector(1,n); */
/* 	ptt=vector(1,n); */
/* 	xit=vector(1,n); */
/* 	*fret=(*func)(p); */
/* 	for (j=1;j<=n;j++) pt[j]=p[j]; */
/* 	for (*iter=1;;++(*iter)) { */
/* 		fp=(*fret); */
/* 		ibig=0; */
/* 		del=0.0; */
/* 		for (i=1;i<=n;i++) { */
/* 			for (j=1;j<=n;j++) xit[j]=xi[j][i]; */
/* 			fptt=(*fret); */
/* 			linmin(p,xit,n,fret,func); */
/* 			if (fabs(fptt-(*fret)) > del) { */
/* 				del=fabs(fptt-(*fret)); */
/* 				ibig=i; */
/* 			} */
/* 		} */
/* 		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) { */
/* 			free_vector(xit,1,n); */
/* 			free_vector(ptt,1,n); */
/* 			free_vector(pt,1,n); */
/* 			return; */
/* 		} */
/* 		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations."); */
/* 		for (j=1;j<=n;j++) { */
/* 			ptt[j]=2.0*p[j]-pt[j]; */
/* 			xit[j]=p[j]-pt[j]; */
/* 			pt[j]=p[j]; */
/* 		} */
/* 		fptt=(*func)(ptt); */
/* 		if (fptt < fp) { */
/* 			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt); */
/* 			if (t < 0.0) { */
/* 				linmin(p,xit,n,fret,func); */
/* 				for (j=1;j<=n;j++) { */
/* 					xi[j][ibig]=xi[j][n]; */
/* 					xi[j][n]=xit[j]; */
/* 				} */
/* 			} */
/* 		} */
/* 	} */
/* } */
/* #undef ITMAX */
/* #undef NRANSI */

/* #define NRANSI */
/* #define TOL 2.0e-4 */

/* int ncom; */
/* double *pcom,*xicom,(*nrfunc)(double []); */

/* void linmin(double p[], double xi[], int n, double *fret, double (*func)(double [])) */
/* { */
/* 	float brent(double ax, double bx, double cx, */
/* 		double (*f)(double), double tol, double *xmin); */
/* 	double f1dim(double x); */
/* 	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, */
/* 		double *fc, double (*func)(double)); */
/* 	int j; */
/* 	double xx,xmin,fx,fb,fa,bx,ax; */

/* 	ncom=n; */
/* 	pcom=vector(1,n); */
/* 	xicom=vector(1,n); */
/* 	nrfunc=func; */
/* 	for (j=1;j<=n;j++) { */
/* 		pcom[j]=p[j]; */
/* 		xicom[j]=xi[j]; */
/* 	} */
/* 	ax=0.0; */
/* 	xx=1.0; */
/* 	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim); */
/* 	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin); */
/* 	for (j=1;j<=n;j++) { */
/* 		xi[j] *= xmin; */
/* 		p[j] += xi[j]; */
/* 	} */
/* 	free_vector(xicom,1,n); */
/* 	free_vector(pcom,1,n); */
/* } */
/* #undef TOL */
/* #undef NRANSI */


/* #define NRANSI */
/* #define ITMAX 100 */
/* #define CGOLD 0.3819660 */
/* #define ZEPS 1.0e-10 */
/* #define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); */

/* double brent(double ax, double bx, double cx, double (*f)(double), double tol, */
/* 	double *xmin) */
/* { */
/* 	int iter; */
/* 	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm; */
/* 	double e=0.0; */

/* 	a=(ax < cx ? ax : cx); */
/* 	b=(ax > cx ? ax : cx); */
/* 	x=w=v=bx; */
/* 	fw=fv=fx=(*f)(x); */
/* 	for (iter=1;iter<=ITMAX;iter++) { */
/* 		xm=0.5*(a+b); */
/* 		tol2=2.0*(tol1=tol*fabs(x)+ZEPS); */
/* 		if (fabs(x-xm) <= (tol2-0.5*(b-a))) { */
/* 			*xmin=x; */
/* 			return fx; */
/* 		} */
/* 		if (fabs(e) > tol1) { */
/* 			r=(x-w)*(fx-fv); */
/* 			q=(x-v)*(fx-fw); */
/* 			p=(x-v)*q-(x-w)*r; */
/* 			q=2.0*(q-r); */
/* 			if (q > 0.0) p = -p; */
/* 			q=fabs(q); */
/* 			etemp=e; */
/* 			e=d; */
/* 			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) */
/* 				d=CGOLD*(e=(x >= xm ? a-x : b-x)); */
/* 			else { */
/* 				d=p/q; */
/* 				u=x+d; */
/* 				if (u-a < tol2 || b-u < tol2) */
/* 					d=SIGN(tol1,xm-x); */
/* 			} */
/* 		} else { */
/* 			d=CGOLD*(e=(x >= xm ? a-x : b-x)); */
/* 		} */
/* 		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d)); */
/* 		fu=(*f)(u); */
/* 		if (fu <= fx) { */
/* 			if (u >= x) a=x; else b=x; */
/* 			SHFT(v,w,x,u) */
/* 			SHFT(fv,fw,fx,fu) */
/* 		} else { */
/* 			if (u < x) a=u; else b=u; */
/* 			if (fu <= fw || w == x) { */
/* 				v=w; */
/* 				w=u; */
/* 				fv=fw; */
/* 				fw=fu; */
/* 			} else if (fu <= fv || v == x || v == w) { */
/* 				v=u; */
/* 				fv=fu; */
/* 			} */
/* 		} */
/* 	} */
/* 	nrerror("Too many iterations in brent"); */
/* 	*xmin=x; */
/* 	return fx; */
/* } */

/* #undef ITMAX */
/* #undef CGOLD */
/* #undef ZEPS */
/* #undef SHFT */
/* #undef NRANSI */



/* #define NRANSI */
/* extern int ncom; */
/* extern double *pcom,*xicom,(*nrfunc)(double []); */

/* double f1dim(double x) */
/* { */
/* 	int j; */
/* 	double f,*xt; */

/* 	xt=vector(1,ncom); */
/* 	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j]; */
/* 	f=(*nrfunc)(xt); */
/* 	free_vector(xt,1,ncom); */
/* 	return f; */
/* } */
/* #undef NRANSI */
