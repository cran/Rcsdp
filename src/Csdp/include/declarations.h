/*
 * CSDPDECLARATIONS is used to prevent redefinitions if this file is included
 * twice.
 */

#ifndef CSDPDECLARATIONS
#define CSDPDECLARATIONS 

/*
  Other important includes that we need.
 */

#include "index.h"
#include "blockmat.h"
#include "parameters.h"

/*
  Our own routines.
  */

void triu(struct blockmatrix A);

void store_packed(struct blockmatrix A, struct blockmatrix B);

void store_unpacked(struct blockmatrix A, struct blockmatrix B);

void alloc_mat_packed(struct blockmatrix A, struct blockmatrix *pB);

void free_mat_packed(struct blockmatrix A);

int structnnz(int n, int k, struct blockmatrix C, struct constraintmatrix *constraints);

int actnnz(int n, int lda, double *A);

int bandwidth(int n, int lda, double *A);

void qreig(int n, double *maindiag, double *offdiag);

void sort_entries(int k, struct blockmatrix C, struct constraintmatrix *constraints);

double norm2(int n, double *x);

double norm1(int n, double *x);

double norminf(int n, double *x);

double Fnorm(struct blockmatrix A);

double Knorm(struct blockmatrix A);

double mat1norm(struct blockmatrix A);

double matinfnorm(struct blockmatrix A);

double calc_pobj(struct blockmatrix C, struct blockmatrix X, 
		 double constant_offset);

double calc_dobj(int k, double *a, double *y, double constant_offset);

double trace_prod(struct blockmatrix A, struct blockmatrix B);

double linesearch(int n, struct blockmatrix dX,
		  struct blockmatrix work1, struct blockmatrix work2, 
		  struct blockmatrix work3, struct blockmatrix cholinv, 
		  double *q, double *z, double *workvec,
		  double stepfrac,double start, int printlevel);

double pinfeas(int k, struct constraintmatrix *constraints,
	       struct blockmatrix X, double *a, double *workvec);

double dinfeas(int k, struct blockmatrix C, 
	       struct constraintmatrix *constraints, double *y, 
	       struct blockmatrix Z, struct blockmatrix work1);

double dimacserr3(int k, struct blockmatrix C, 
	       struct constraintmatrix *constraints, double *y, 
	       struct blockmatrix Z, struct blockmatrix work1);

void op_a(int k, struct constraintmatrix *constraints,
	  struct blockmatrix X, double *result);

void op_at(int k, double *y, struct constraintmatrix *constraints,
	   struct blockmatrix result);

void makefill(int k, struct blockmatrix C, 
	      struct constraintmatrix *constraints, 
	      struct constraintmatrix *pfill, struct blockmatrix work1, 
	      int printlevel);

void op_o(int k, struct constraintmatrix *constraints,
	  struct sparseblock **byblocks, struct blockmatrix Zi, 
          struct blockmatrix X, double *O, struct blockmatrix work1, 
          struct blockmatrix work2);

void addscaledmat(struct blockmatrix A, double scale, struct blockmatrix B,
		  struct blockmatrix C);

void zero_mat(struct blockmatrix A);

void add_mat(struct blockmatrix A,struct blockmatrix B);

void sym_mat(struct blockmatrix A);

void make_i(struct blockmatrix A);

void copy_mat(struct blockmatrix A, struct blockmatrix B);

void mat_mult(double scale1, double scale2, struct blockmatrix A,
	      struct blockmatrix B, struct blockmatrix C);

void mat_multspa(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);

void mat_multspb(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);

void mat_multspc(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);

void mat_mult_raw(int n, double scale1, double scale2, double *ap,
		  double *bp, double *cp);

void mat_mult_rawatlas(int n, double scale1, double scale2, double *ap,
		  double *bp, double *cp);

void matvec(struct blockmatrix A, double *x, double *y);

void alloc_mat(struct blockmatrix A, struct blockmatrix *pB);

void free_mat(struct blockmatrix A);

void initparams(struct paramstruc *params, int *pprintlevel);

void initsoln(int n, int k, struct blockmatrix C, double *a, 
	      struct constraintmatrix *constraints, struct blockmatrix *pX0,
	      double **py0, struct blockmatrix *pZ0);

void trans(struct blockmatrix A);

void chol_inv(struct blockmatrix A, struct blockmatrix B);

int chol(struct blockmatrix A);

int solvesys(int m, int ldam, double *A, double *rhs);

int user_exit(int n, int k, struct blockmatrix C, double *a, double dobj,
	      double pobj, double constant_offset, 
	      struct constraintmatrix *constraints, struct blockmatrix X,
	      double *y, struct blockmatrix Z, struct paramstruc params);

int read_sol(char *fname, int n, int k, struct blockmatrix C, 
	     struct blockmatrix *pX, double **py, struct blockmatrix *pZ);

int read_prob(char *fname, int *pn, int *pk, struct blockmatrix *pC,
	      double **pa, struct constraintmatrix **pconstraints,
	      int printlevel);

int write_prob(char *fname, int n, int k, struct blockmatrix C,
	       double *a, struct constraintmatrix *constraints);

int write_sol(char *fname, int n, int k, struct blockmatrix X,
	      double *y, struct blockmatrix Z);

void free_prob(int n, int k, struct blockmatrix C, double *a, 
	       struct constraintmatrix *constraints, struct blockmatrix X,
	       double *y, struct blockmatrix Z);

int sdp(int n, int k, struct blockmatrix C, double *a, double constant_offset,
	struct constraintmatrix *constraints, struct sparseblock **byblocks,
	struct constraintmatrix fill, struct blockmatrix X, double *y, 
	struct blockmatrix Z, struct blockmatrix cholxinv,
	struct blockmatrix cholzinv, double *pobj, double *dobj, 
	struct blockmatrix work1, struct blockmatrix work2,     
	struct blockmatrix work3,
	double *workvec1, double *workvec2,
	double *workvec3, double *workvec4, double *workvec5,
	double *workvec6, double *workvec7, double *workvec8, 
	double *diagO, struct blockmatrix bestx, double *besty,
	struct blockmatrix bestz, struct blockmatrix Zi, double *O,
	double *rhs, struct blockmatrix dZ, struct blockmatrix dX, 
	double *dy, double *dy1, double *Fp, int printlevel, 
	struct paramstruc parameters);

int easy_sdp(int n, int k, struct blockmatrix C, double *a, 
	     struct constraintmatrix *constraints, double constant_offset,
	     struct blockmatrix *pX, double **py, struct blockmatrix *pZ,
	     double *ppobj, double *pdobj);

void tweakgap(int n, int k, double *a, struct constraintmatrix *constraints,
	      double gap, struct blockmatrix Z, struct blockmatrix dZ, 
	      double *y, double *dy, struct blockmatrix work1, 
	      struct blockmatrix work2, struct blockmatrix work3, 
	      struct blockmatrix work4, double *workvec1, double *workvec2,
	      double *workvec3, double *workvec4, int printlevel);

int bisect_(int *n, double *eps1, double *d, double *e, double *e2,
	    double *lb, double *ub, int *mm, int *m, double *w, int *ind, 
	    int *ierr, double *rv4, double *rv5);

/*
  BLAS and LINPACK stuff.
  */

/*
  First, BLAS routines.
 */

#ifdef CAPSBLAS
#ifdef NOUNDERBLAS
double DNRM2(int *, double *, int *);
double DASUM(int *, double *, int *);
double DDOT(int *, double *, int *, double *, int *);
int IDAMAX(int *, double*, int *);
void DGEMM(const char *, const char *, int *,   int *, int *, double *,   double *, int *, double *, int *, double *, double *,  int *);
void DGEMV(const char *, int *, int *,  double *, double *, int *,  double *, int *, double *,  double*, int *);
#else
double DNRM2_(int *, double *, int *);
double DASUM_(int *, double *, int *);
double DDOT_(int *, double *, int *, double *, int *);
int IDAMAX_(int *, double*, int *);
void DGEMM_(const char *, const char *, int *,   int *, int *, double *,   double *, int *, double *, int *, double *, double *,  int *);
void DGEMV_(const char *, int *, int *,  double *, double *, int *,   *, int *, double *,  double*, int *);
#endif
#else
#ifdef NOUNDERBLAS
double dnrm2(int *, double *, int *);
double dasum(int *, double *, int *);
double ddot(int *, double *, int *, double *, int *);
int idamax(int *, double*, int *);
void dgemm(const char *, const char *, int *,   int *, int *, double *,   double *, int *, double *, int *, double *, double *,  int *);
void dgemv(const char *, int *, int *,  double *, double *, int *,  double *, int *, double *,  double*, int *);
#else
double dnrm2_(int *, double *, int *);
double dasum_(int *, double *, int *);
double ddot_(int *, double *, int *, double *, int *);
int idamax_(int *, double*, int *);
void dgemm_(const char *, const char *, int *,   int *, int *, double *,   double *, int *, double *, int *, double *, double *,  int *);
void dgemv_(const char *, int *, int *,  double *, double *, int *,  double *, int *, double *,  double*, int *);
#endif
#endif

/*
  LAPACK next.
 */

#ifdef CAPSLAPACK
#ifdef NOUNDERLAPACK
void DPOTRF(const char *, const int *, double *, const int *, int *);
void DPOTRS(const char *, int *, int *, double *, int *, double *, int *, int *);
void DTRTRI(const char *, const char *, int *, double *, int *, int *);
#else
void DPOTRF_(const char *, const int *, double *, const int *, int *);
void DPOTRS_(const char *, int *, int *, double *, int *, double *, int *, int *);
void DTRTRI_(const char *, const char *, int *, double *, int *, int *);
#endif
#else
#ifdef NOUNDERLAPACK
void dpotrf(const char *, const int *, double *, const int *, int *);
void dpotrs(const char *, int *, int *, double *, int *, double *, int *, int *);
void dtrtri(const char *, const char *, int *, double *, int *, int *);
#else
void dpotrf_(const char *, const int *, double *, const int *, int *);
void dpotrs_(const char *, int *, int *, double *, int *, double *, int *, int *);
void dtrtri_(const char *, const char *, int *, double *, int *, int *);
#endif
#endif


#endif
