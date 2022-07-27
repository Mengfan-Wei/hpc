static char help[] = "Illustrate the use of MatResetPreallocation.\n";

#include "petscmat.h"

int main(int argc,char **argv)
{
  Mat             A;
  MPI_Comm        comm;
  PetscMPIInt     rank;
  PetscInt        n=5,m=5,*dnnz,*onnz,i,rstart,rend,M,N;
  PetscErrorCode  ierr;

  ierr = PetscInitialize(&argc,&argv,0,help);if (ierr) return ierr;
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  ierr = PetscMalloc2(m,&dnnz,m,&onnz);CHKERRQ(ierr);
  for (i=0; i<m; i++) 
  {
    dnnz[i] = 1;
    onnz[i] = 1;
  }
  ierr = MatCreateAIJ(comm,m,n,PETSC_DETERMINE,PETSC_DETERMINE,PETSC_DECIDE,dnnz,PETSC_DECIDE,onnz,&A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = PetscFree2(dnnz,onnz);CHKERRQ(ierr);

  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&rstart,&rend);CHKERRQ(ierr);
  ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "Matrix size is %d by %d \n", M, N);
  ierr = PetscPrintf(PETSC_COMM_SELF, "rank [%d] rstart = %d , rend = %d \n", rank, rstart, rend);
  for (i=rstart; i<rend; i++) 
  {
    ierr = MatSetValue(A,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
    if (rend<N) 
      ierr = MatSetValue(A,i,rend,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

// EOF
