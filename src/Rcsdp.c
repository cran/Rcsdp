/*
 * R interface to CSDP semidefinite programming library
 *
 * Created: Hector Corrada Bravo <hcorrada@gmail.com>
 * February 22, 2008
 */

#include <R.h>
#include <Rinternals.h>
#include <declarations.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP int_vector_csdp2R(int, int*);
SEXP double_vector_csdp2R(int, double*);
int *int_vector_R2csdp(int, SEXP);
double *double_vector_R2csdp(int, SEXP);
struct blockmatrix blkmatrix_R2csdp(SEXP);
SEXP blkmatrix_csdp2R(struct blockmatrix);
SEXP constraints_csdp2R(int, struct constraintmatrix *);
struct constraintmatrix *constraints_R2csdp(SEXP);

void free_constraints(int, struct constraintmatrix *);

SEXP csdp(SEXP n_p,
	  SEXP nconstraints_p,
	  SEXP nblocks_p,
	  SEXP blocktypes_p,
	  SEXP blocksizes_p,
	  SEXP C_p,
	  SEXP A_p,
	  SEXP b_p)
{
  SEXP X_p, Z_p, y_p, ret, pobj_p, dobj_p, status_p;
  enum AIJ_SLOTS {AIJ_NNZ, AIJ_IIND, AIJ_JIND, AIJ_X};

  int n, nblocks, nconstraints, *blocktypes, *blocksizes;

  struct blockmatrix C;
  struct constraintmatrix *constraints;
  struct blockmatrix X, Z;
  double *y, *b;
  double pobj, dobj;
  int status;

  n = INTEGER(n_p)[0];
  nblocks = INTEGER(nblocks_p)[0];
  nconstraints = INTEGER(nconstraints_p)[0];
  blocktypes = INTEGER(blocktypes_p);
  blocksizes = INTEGER(blocksizes_p);

  /*
   * setup C
   */
  C = blkmatrix_R2csdp(C_p);

  /*
   * setup constraints 
   */
  constraints = constraints_R2csdp(A_p);

  /*
   * Allocate storage for RHS
   */
  b = double_vector_R2csdp(nconstraints,b_p);
  if (b==NULL) error("Failed to allocate storage for RHS vector b!\n");

  /*
   * Create an initial solution. This allocates space for X, y, and Z,
   * and sets initial values
   */
  initsoln(n,nconstraints,C,b,constraints,&X,&y,&Z);

  /*
   * Solve the problem
   */
  status = easy_sdp(n,nconstraints,C,b,constraints,0.0,&X,&y,&Z,&pobj,&dobj);

  /* 
   * Grab the results
   */
  
  /*
   * Grab X
   */
  X_p = PROTECT(blkmatrix_csdp2R(X));
  
  /*
   * Grab Z
   */
  Z_p = PROTECT(blkmatrix_csdp2R(Z));
  
  /* Copy y */
  y_p = PROTECT(double_vector_csdp2R(nconstraints, y));

  pobj_p = PROTECT(allocVector(REALSXP,1)); REAL(pobj_p)[0] = pobj;
  dobj_p = PROTECT(allocVector(REALSXP,1)); REAL(dobj_p)[0] = dobj;
  status_p = PROTECT(allocVector(INTSXP,1)); INTEGER(status_p)[0] = status;
  
  free_prob(n,nconstraints,C,b,constraints,X,y,Z);

  ret = PROTECT(allocVector(VECSXP,6));
  SET_VECTOR_ELT(ret,0,X_p);
  SET_VECTOR_ELT(ret,1,Z_p);
  SET_VECTOR_ELT(ret,2,y_p);
  SET_VECTOR_ELT(ret,3,pobj_p);
  SET_VECTOR_ELT(ret,4,dobj_p);
  SET_VECTOR_ELT(ret,5,status_p);

  UNPROTECT(7);
  return ret;
}

SEXP get_prob_info(struct blockmatrix X)
{
  SEXP ret, types, sizes;
  int nblocks, i;
  int *intvec;

  nblocks = X.nblocks;
  PROTECT(ret = allocVector(VECSXP,2));

  PROTECT(types = allocVector(INTSXP,nblocks+1));
  intvec = INTEGER(types);
  for (i=1; i<=nblocks; i++)
    intvec[i] = X.blocks[i].blockcategory == MATRIX ? 1 : 2;
  SET_VECTOR_ELT(ret, 0, types);
  
  PROTECT(sizes = allocVector(INTSXP,nblocks+1));
  intvec = INTEGER(sizes);
  for (i=1; i<=nblocks; i++)
    intvec[i] = X.blocks[i].blocksize;
  SET_VECTOR_ELT(ret, 1, sizes);

  UNPROTECT(3);
  return ret;
}

SEXP readsdpa(SEXP filename,
	      SEXP verbose)
{
  int n,k;
  struct blockmatrix C;
  double *b;
  struct constraintmatrix *constraints;

  int printlevel;
  int status;
  char *fname;
  SEXP ret;

  fname = (char *) CHAR(STRING_ELT(filename,0));
  printlevel = INTEGER(verbose)[0];
  status = read_prob(fname,&n,&k,&C,&b,&constraints,printlevel);
  if (status) error("Error reading sdpa file %s, status:%d\n",fname,status);

  ret = PROTECT(allocVector(VECSXP,4));
  SET_VECTOR_ELT(ret,0,blkmatrix_csdp2R(C));
  SET_VECTOR_ELT(ret,1,constraints_csdp2R(k,constraints));
  SET_VECTOR_ELT(ret,2,double_vector_csdp2R(k,b));
  SET_VECTOR_ELT(ret,3,get_prob_info(C));

  free(b);
  free_mat(C);
  free_constraints(k,constraints);

  UNPROTECT(1);
  return ret;
}

SEXP readsdpa_sol(SEXP filename,
		  SEXP n_p,
		  SEXP k_p,
		  SEXP C_p)
{
  char *fname;
  int n, k, status;
  struct blockmatrix C;
  struct blockmatrix X;
  struct blockmatrix Z;
  double *y;

  SEXP ret, X_p, y_p, Z_p;

  n = INTEGER(n_p)[0];
  k = INTEGER(k_p)[0];
  C = blkmatrix_R2csdp(C_p);

  fname = (char *) CHAR(STRING_ELT(filename,0));
  status = read_sol(fname,n,k,C,&X,&y,&Z);
  if (status) {
    free_mat(C);
    free_mat(X);
    free(y);
    free_mat(Z);
    error("Reading reading solution in file %s: status %d\n",fname,status);
  }

  ret = PROTECT(allocVector(VECSXP,3));
  X_p = PROTECT(blkmatrix_csdp2R(X));
  y_p = PROTECT(double_vector_csdp2R(k,y));
  Z_p = PROTECT(blkmatrix_csdp2R(Z));

  free_mat(C);
  free_mat(X);
  free(y);
  free_mat(Z);

  SET_VECTOR_ELT(ret, 0, X_p);
  SET_VECTOR_ELT(ret, 1, y_p);
  SET_VECTOR_ELT(ret, 2, Z_p);
  UNPROTECT(4);
  return ret;
}
		  
SEXP writesdpa(SEXP filename,
	       SEXP n_p,
	       SEXP nconstraints_p,
	       SEXP nblocks_p,
	       SEXP blocktypes_p,
	       SEXP blocksizes_p,
	       SEXP C_p,
	       SEXP A_p,
	       SEXP b_p)
{
  struct blockmatrix C;
  struct constraintmatrix *constraints;
  double *b;

  int n, nblocks, nconstraints, status;
  char *fname;

  SEXP ret;

  n = INTEGER(n_p)[0];
  nblocks = INTEGER(nblocks_p)[0];
  nconstraints = INTEGER(nconstraints_p)[0];

  fname = (char *) CHAR(STRING_ELT(filename,0));  

  /*
   * setup C
   */
  C = blkmatrix_R2csdp(C_p);

  /*
   * setup constraints 
   */
  constraints = constraints_R2csdp(A_p);

  /*
   * Allocate storage for RHS
   */
  b = double_vector_R2csdp(nconstraints,b_p);
  if (b==NULL) error("Failed to allocate storage for RHS vector b!\n");

  status = write_prob(fname,n,nconstraints,C,b,constraints);
  free_mat(C);
  free_constraints(nconstraints, constraints);
  free(b);

  PROTECT(ret = allocVector(INTSXP,1));
  INTEGER(ret)[0] = status;
  UNPROTECT(1);
  return ret;
}

SEXP writesdpa_sol(SEXP filename,
		   SEXP n_p,
		   SEXP k_p,
		   SEXP X_p,
		   SEXP y_p,
		   SEXP Z_p)
{
  int n,k,status;
  struct blockmatrix X;
  struct blockmatrix Z;
  double *y;
  char *fname;

  SEXP ret;

  n = INTEGER(n_p)[0];
  k = INTEGER(k_p)[0];

  fname = (char *) CHAR(STRING_ELT(filename,0));
  X = blkmatrix_R2csdp(X_p);
  Z = blkmatrix_R2csdp(Z_p);
  y = double_vector_R2csdp(k,y_p);

  status = write_sol(fname,n,k,X,y,Z);
  free_mat(X);
  free(y);
  free_mat(Z);

  PROTECT(ret = allocVector(INTSXP,1));
  INTEGER(ret)[0] = status;
  UNPROTECT(1);
  return ret;

}


void free_constraints(int k,
		      struct constraintmatrix *constraints)
{
  int i;
  struct sparseblock *ptr;
  struct sparseblock *oldptr;

  
    if (constraints != NULL)
    {
      for (i=1; i<=k; i++)
	{
	  /*
	   * Get rid of constraint i.
	   */
	  
	  ptr=constraints[i].blocks;
	  while (ptr != NULL)
	    {
	      free(ptr->entries);
	      free(ptr->iindices);
	      free(ptr->jindices);
	      oldptr=ptr;
	      ptr=ptr->next;
	      free(oldptr);
	    };
	};
      /*
       * Finally, free the constraints array.
       */

      free(constraints);
    };
}


static const R_CallMethodDef CallEntries[] = {
  {"csdp",          (DL_FUNC) &csdp,          8},
  {"readsdpa",      (DL_FUNC) &readsdpa,      2},
  {"readsdpa_sol",  (DL_FUNC) &readsdpa_sol,  4},
  {"writesdpa",     (DL_FUNC) &writesdpa,     9},
  {"writesdpa_sol", (DL_FUNC) &writesdpa_sol, 6},
  {NULL, NULL, 0}
};

void R_init_Rcsdp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
