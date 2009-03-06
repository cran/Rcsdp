#include <R.h>
#include <Rinternals.h>
#include <declarations.h>

SEXP int_vector_csdp2R(int n,
		       int *y)
{
  SEXP ret;
  int i;
  int *intvec;

  PROTECT(ret = allocVector(INTSXP,n+1));
  intvec = INTEGER(ret);
  for (i=1; i<=n; i++)
    intvec[i] = y[i];
  return ret;
}

SEXP double_vector_csdp2R(int n,
			  double *y)
{
  SEXP ret;
  int i;
  double *dblvec;

  PROTECT(ret = allocVector(REALSXP,n+1));
  dblvec = REAL(ret);
  for (i=1; i<=n; i++)
    dblvec[i] = y[i];
  return ret;
}

int * int_vector_R2csdp(int n,
			SEXP y)
{
  int *ret;
  int *intvec;
  int i;

  ret = (int *) malloc((n+1) * sizeof(int));
  if (ret == NULL)
    return NULL;

  intvec = INTEGER(y);
  for (i=1; i<=n; i++)
    ret[i] = intvec[i];
  return ret;
}

double * double_vector_R2csdp(int n,
			      SEXP y)
{
  double *ret;
  double *dblvec;
  int i;

  ret = (double *) malloc((n+1) * sizeof(double));
  if (ret == NULL)
    return NULL;

  dblvec = REAL(y);
  for (i=1; i<=n; i++)
    ret[i] = dblvec[i];
  return ret;
}

struct blockmatrix blkmatrix_R2csdp(SEXP X)
{
  struct blockmatrix ret;
  SEXP blocks, cur_block;
  int nblocks, blksize, blktype, allocsize;
  int j, k;
  double *dblvec;

  enum blkmatrix_slots {NBLOCKS, BLOCKS};
  enum blkrec_slots {BLOCKSIZE, BLOCKCATEGORY,DATA};

  nblocks = INTEGER(VECTOR_ELT(X,NBLOCKS))[0];
  blocks = VECTOR_ELT(X,BLOCKS);

  ret.nblocks = nblocks;
  ret.blocks = (struct blockrec *) malloc((nblocks + 1) * sizeof(struct blockrec));
  if (ret.blocks == NULL)
    error("Error allocating blkmatrix blocks");

  for (j=1; j<=nblocks; j++) {
    cur_block = VECTOR_ELT(blocks,j-1);
    blksize = INTEGER(VECTOR_ELT(cur_block,BLOCKSIZE))[0];
    ret.blocks[j].blocksize = blksize;

    blktype = INTEGER(VECTOR_ELT(cur_block,BLOCKCATEGORY))[0];
    ret.blocks[j].blockcategory = (blktype == 1) ? MATRIX : DIAG;
    
    if (blktype == MATRIX) {
      allocsize = blksize*blksize;
      ret.blocks[j].data.mat = (double *) malloc(allocsize * sizeof(double));
      if (ret.blocks[j].data.mat == NULL)
	error("Error allocating block matrix data, s block");

      dblvec = REAL(VECTOR_ELT(cur_block,DATA));
      for (k=0; k<allocsize; k++)
	ret.blocks[j].data.mat[k] = dblvec[k];
    }
    else {
      ret.blocks[j].data.vec = double_vector_R2csdp(blksize, VECTOR_ELT(cur_block,DATA));
      if (ret.blocks[j].data.vec == NULL)
	error("Error allocating block matrix data, l block");
    }
  }
  return ret;
}

SEXP blkmatrix_csdp2R(struct blockmatrix X)
{
  SEXP ret;
  SEXP blocks, nblocks, cur_block;
  SEXP blocksize, blockcategory, data;

  int j,k, allocsize;
  double *dblvec;

  enum blkmatrix_slots {NBLOCKS, BLOCKS};
  enum blkrec_slots {BLOCKSIZE, BLOCKCATEGORY,DATA};

  PROTECT(ret = allocVector(VECSXP, 2));

  PROTECT(nblocks = allocVector(INTSXP, 1));
  INTEGER(nblocks)[0] = X.nblocks;
  SET_VECTOR_ELT(ret, NBLOCKS, nblocks);

  PROTECT(blocks = allocVector(VECSXP, X.nblocks));
  for (j=1; j<=X.nblocks; j++) {
    PROTECT(cur_block = allocVector(VECSXP,3));
    PROTECT(blocksize = allocVector(INTSXP,1));
    INTEGER(blocksize)[0] = X.blocks[j].blocksize;
    
    PROTECT(blockcategory = allocVector(INTSXP,1));
    INTEGER(blockcategory)[0] = (X.blocks[j].blockcategory == MATRIX) ? 1 : 2;

    if (X.blocks[j].blockcategory == MATRIX) {
      allocsize = X.blocks[j].blocksize * X.blocks[j].blocksize;
      PROTECT(data = allocVector(REALSXP,allocsize));
      dblvec = REAL(data);
      for (k=0; k<allocsize; k++)
	dblvec[k] = X.blocks[j].data.mat[k];
    }
    else
      data = double_vector_csdp2R(X.blocks[j].blocksize, X.blocks[j].data.vec);

    SET_VECTOR_ELT(cur_block, BLOCKSIZE, blocksize);
    SET_VECTOR_ELT(cur_block, BLOCKCATEGORY, blockcategory);
    SET_VECTOR_ELT(cur_block, DATA, data);
    SET_VECTOR_ELT(blocks,j-1,cur_block);
    UNPROTECT(4);
  }
  SET_VECTOR_ELT(ret, BLOCKS, blocks);
  UNPROTECT(2);
  return ret;
}

SEXP constraints_csdp2R(int numconstraints,
		      struct constraintmatrix *constraints)
{
  struct sparseblock *ptr;
  int i, j, k;
  SEXP ret, Ai, Aij;
  SEXP iindices, jindices, entries;
  SEXP blocknum, blocksize, constraintnum;
  SEXP numentries;

  enum constraint_blockslots {IIND, JIND, ENTRIES, BLOCKNUM, BLOCKSIZE, CONSTRAINTNUM, NUMENTRIES};
 
  int *intvec;
  double *dblvec;
  int nblocks, nnz;

  PROTECT(ret = allocVector(VECSXP, numconstraints));

  if (constraints != NULL) {
    for (i=1; i<=numconstraints; i++) {
      ptr = constraints[i].blocks;
      nblocks = 0;
      while (ptr != NULL) {
	nblocks += 1;
	ptr = ptr->next;
      }
      PROTECT(Ai = allocVector(VECSXP,nblocks));

      ptr = constraints[i].blocks;
      for (j=1; j<=nblocks; j++) {
	PROTECT(Aij = allocVector(VECSXP,7));
	nnz = ptr->numentries;

	PROTECT(numentries = allocVector(INTSXP,1));
	INTEGER(numentries)[0] = nnz;
	SET_VECTOR_ELT(Aij, NUMENTRIES, numentries);

	PROTECT(blocknum = allocVector(INTSXP,1));
	INTEGER(blocknum)[0] = ptr->blocknum;
	SET_VECTOR_ELT(Aij, BLOCKNUM, blocknum);

	PROTECT(blocksize = allocVector(INTSXP,1));
	INTEGER(blocksize)[0] = ptr->blocksize;
	SET_VECTOR_ELT(Aij, BLOCKSIZE, blocksize);
	
	PROTECT(constraintnum = allocVector(INTSXP,1));
	INTEGER(constraintnum)[0] = ptr->constraintnum;
	SET_VECTOR_ELT(Aij, CONSTRAINTNUM, constraintnum);

	iindices = int_vector_csdp2R(nnz, ptr->iindices);
	SET_VECTOR_ELT(Aij,IIND, iindices);

	jindices = int_vector_csdp2R(nnz, ptr->jindices);
	SET_VECTOR_ELT(Aij, JIND, jindices);

	entries = double_vector_csdp2R(nnz, ptr->entries);
	SET_VECTOR_ELT(Aij, ENTRIES, entries);

	SET_VECTOR_ELT(Ai, j-1, Aij);
	UNPROTECT(8);
	ptr = ptr->next;
      }
      SET_VECTOR_ELT(ret, i-1, Ai);
      UNPROTECT(1);
    }
  }
  return ret;
}

struct constraintmatrix *constraints_R2csdp(SEXP A)
{
  struct constraintmatrix *constraints;
  struct sparseblock *blockptr;

  int nconstraints, nblocks;
  SEXP Ai, Aij;

  int i,j;

  enum constraint_blockslots {IINDICES, JINDICES, ENTRIES, BLOCKNUM, BLOCKSIZE, CONSTRAINTNUM, NUMENTRIES};

  nconstraints = LENGTH(A);
  constraints = (struct constraintmatrix *) malloc((nconstraints + 1) * sizeof(struct constraintmatrix));
  if (constraints == NULL) error("Failed to allocate storage for constraints!\n");
    
  for (i=1; i<=nconstraints; i++) {
    Ai = VECTOR_ELT(A,i-1);
    /*
     * Terminate block linked list with NULL
     */
    constraints[i].blocks = NULL;
    nblocks = LENGTH(Ai);

    for (j=nblocks; j>=1; j--) {
      Aij = VECTOR_ELT(Ai,j-1);
      
      /*
       * Allocate block data structure
       */
      blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));
      if (blockptr == NULL) error("Allocation of constraint block failed!\n");
      
      /*
       * Initialize block data structure
       */
      blockptr->blocknum=INTEGER(VECTOR_ELT(Aij,BLOCKNUM))[0];
      blockptr->blocksize=INTEGER(VECTOR_ELT(Aij,BLOCKSIZE))[0];
      blockptr->constraintnum=INTEGER(VECTOR_ELT(Aij,CONSTRAINTNUM))[0];
      blockptr->next=NULL;
      blockptr->nextbyblock=NULL;
      blockptr->numentries=INTEGER(VECTOR_ELT(Aij,NUMENTRIES))[0];

      /*
       * Enter data
       */
      blockptr->iindices = int_vector_R2csdp(blockptr->numentries, VECTOR_ELT(Aij,IINDICES));
      if (blockptr->iindices == NULL) error("Allocation of constraint block failed\n");

      blockptr->jindices = int_vector_R2csdp(blockptr->numentries, VECTOR_ELT(Aij,JINDICES));
      if (blockptr->jindices == NULL) error("Allocation of constraint block failed\n");

      blockptr->entries = double_vector_R2csdp(blockptr->numentries, VECTOR_ELT(Aij,ENTRIES));
      if (blockptr->entries == NULL) error("Allocation of constraint block failed\n");

      /*
       * Insert block into linked list
       */
      blockptr->next=constraints[i].blocks;
      constraints[i].blocks=blockptr;
    }
  }
  return constraints;
}
