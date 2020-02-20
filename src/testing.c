#include <R.h>
#include <Rinternals.h>
#include <declarations.h>

SEXP int_vector_csdp2R(int, int*);
SEXP double_vector_csdp2R(int, double*);
int *int_vector_R2csdp(int, SEXP);
double *double_vector_R2csdp(int, SEXP);
struct blockmatrix blkmatrix_R2csdp(SEXP);
SEXP blkmatrix_csdp2R(struct blockmatrix);
SEXP constraints_csdp2R(int, struct constraintmatrix *);
struct constraintmatrix *constraints_R2csdp(SEXP);
void free_constraints(int, struct constraintmatrix *);

SEXP test_int_vector(SEXP n_p,
		     SEXP v)
{
  int *vv, n;
  SEXP ret;

  n = INTEGER(n_p)[0];
  vv = int_vector_R2csdp(n,v);
  ret = int_vector_csdp2R(n,vv);
  free(vv);
  return ret;
}

SEXP test_double_vector(SEXP n_p,
			SEXP v)
{
  int n;
  double *vv;
  SEXP ret;

  n = INTEGER(n_p)[0];
  vv = double_vector_R2csdp(n,v);
  ret = double_vector_csdp2R(n,vv);
  free(vv);
  return ret;
  
}

SEXP test_blkmatrix(SEXP X)
{
  struct blockmatrix XX;
  SEXP ret;

  XX = blkmatrix_R2csdp(X);
  ret = blkmatrix_csdp2R(XX);
  free_mat(XX);
  return ret;
}

SEXP test_constraints(SEXP k_p,
		      SEXP A)
{
  int k;
  struct constraintmatrix *AA;
  SEXP ret;

  k = INTEGER(k_p)[0];
  AA = constraints_R2csdp(A);
  ret = constraints_csdp2R(k,AA);
  free_constraints(k,AA);
  return ret;
}
