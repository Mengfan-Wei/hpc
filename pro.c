static char help[] = "Solves a tridiagonal linear system.\n\n";

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

  Note:  The corresponding uniprocessor example is ex1.c
*/
#include <petscksp.h>
int main(int argc,char **args)
{
  Vec            us,u,f;          /* approx solution, RHS, exact solution */
  Mat            A;                /* linear system matrix */
  //KSP            ksp;              /* linear solver context */
  //PC             pc;               /* preconditioner context */
  PetscReal      dx=0.01, dt=0.001, t, L=1.0, k=1.0, r=k*dt/(dx*dx), tol=1000.*PETSC_MACHINE_EPSILON ;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i,n = 10,col[3],rank,rstart,rend,nlocal,max_ite=10000,ite=0;
  PetscScalar    one = 1.0, zero=0.0, value[3],lambda, ;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Az = y.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed. For this simple case let PETSc decide how
     many elements of the vector are stored on each processor. The second
     argument to VecSetSizes() below causes PETSc to decide.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);  
  ierr = VecSetSizes(u,PETSC_DECIDE,n);CHKERRQ(ierr); 
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&us);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&f);CHKERRQ(ierr);

  /* Identify the starting and ending mesh points on each
     processor for the interior part of the mesh. We let PETSc decide
     above. */
  ierr = VecGetOwnershipRange(u,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);
  
  /*
     Set the vector elements.
  */
  ierr = VecSet(z,zero);CHKERRQ(ierr);
  if( rank == 0 )
  { 
 	for i 
  	ierr = VecSetValues(z,1,&i,exp(i),INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);
  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.

     We pass in nlocal as the "local" size of the matrix to force it
     to have the same parallel layout as the vector created above.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /*
     Assemble matrix.

     The linear system is distributed across the processors by
     chunks of contiguous rows, which correspond to contiguous
     sections of the mesh on which the problem is discretized.
     For matrix assembly, each processor contributes entries for
     the part that it owns locally.
  */
  if (!rstart) {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1-2r; value[1] =r;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  if (rend == n) {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = r; value[1] = 1-2r;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Set entries corresponding to the mesh interior */
  value[0] = r; value[1] = 1-2r; value[2] = r;
  for (i=rstart; i<rend; i++) {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


  while ((PetscAbsReal(norm_lag - norm) > tol) && (ite < max_ite)){
	norm = norm_lag;
  	ierr = MatMultAdd(A,u,f,us);CHKERRQ(ierr);
  	ierr = VecNorm(y,NORM_2,&norm_lag);CHKERRQ(ierr);
	ierr = VecScale(y,(PetscScalar)one/norm_lag);CHKERRQ(ierr);
  	ierr = VecCopy(y,z);CHKERRQ(ierr);
        
	if (ite % 20 == 0){
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"Error %g, Iterations %D\n",(double)(norm_lag - norm),ite);CHKERRQ(ierr);
	}
	ite++;
  }

  ierr = MatMultTranspose(A,z,r);CHKERRQ(ierr);
  ierr = VecDot(r,z,&lambda);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"The lagest eigenvaule is %g\n The corresponding eigenvector is\n",(double)lambda);CHKERRQ(ierr);
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  ierr = VecDestroy(&r);CHKERRQ(ierr); ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr); 
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  ierr = PetscFinalize();
  return ierr;
}

// EOF`
