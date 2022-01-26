#include "petscviewer.h"

#include <slepceps.h>
#include <petscsys.h>
#include <petscmat.h>
#include <mpi.h>

void early_exit(int status, const char msg[]) {
    PetscPrintf(PETSC_COMM_WORLD, msg);
    SlepcFinalize();
    exit(status);
}

int main(int argc, char* argv[]) {
    PetscErrorCode          ierr;
    Mat                     A, B;
    Vec                     xr, xi;
    PetscViewer             viewer;
    EPS                     eps;
    ST                      st;
    KSP                     ksp;
    PC                      pc;
    PetscInt                nconv;

    double eigenvalue_r, eigenvalue_i;

    double target_eval = 1.47;

    // ------------------------------------------------------------------------------------------------------------------
    // Basic Setup
    // ------------------------------------------------------------------------------------------------------------------
    int num_threads;
    ierr = SlepcInitialize(&argc, &argv, (char *) 0, NULL); 
    if (ierr) early_exit(1, "Failed to initialize Slepc!");

    MPI_Comm_size(PETSC_COMM_WORLD, &num_threads);
    if (num_threads != 1) early_exit(2, "GEP Solver is only setup to use 1 MPI Thread!");

    // ------------------------------------------------------------------------------------------------------------------
    // Load Matrices From Files 
    // ------------------------------------------------------------------------------------------------------------------
    ierr = MatCreate(PETSC_COMM_WORLD, &A);                         if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = MatSetType(A, MATSEQAIJ);                                if (ierr) early_exit(3, "Unable to Create A Matrix");

    ierr =  MatCreate(PETSC_COMM_WORLD, &B);                        if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = MatSetType(B, MATSEQAIJ);                                if (ierr) early_exit(3, "Unable to Create B Matrix");

    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);            if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = PetscViewerSetType(viewer, "binary");                    if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);          if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = PetscViewerFileSetName(viewer, "./test_a.dat");          if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = MatLoad(A, viewer);                                      if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(3, "Unable to Create A Matrix");

    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);            if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerSetType(viewer, "binary");                    if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);          if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerFileSetName(viewer, "./test_b.dat");          if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = MatLoad(B, viewer);                                      if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(3, "Unable to Create B Matrix");

    ierr = MatCreateVecs(A, NULL, &xr);                             if (ierr) early_exit(3, "Unable to Create Eigenvector");
    ierr = MatCreateVecs(A, NULL, &xi);                             if (ierr) early_exit(3, "Unable to Create Eigenvector");

    // ------------------------------------------------------------------------------------------------------------------
    // Setup Eigenvalue Problem
    // ------------------------------------------------------------------------------------------------------------------
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);                   if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    ierr = EPSSetOperators(eps, A, B);                          if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    ierr = EPSSetProblemType(eps, EPS_GHEP);                    if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    ierr = EPSSetFromOptions(eps);                              if (ierr) early_exit(4, "Unable to Setup Eigenproblem");

    ierr = EPSSetTolerances(eps, 1.0e-15, 100);                 if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    ierr = EPSSetType(eps, "krylovschur");                      if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);    if (ierr) early_exit(4, "Unable to Setup Eigenproblem");
    
    ierr = EPSSetTarget(eps, target_eval);                      if (ierr) early_exit(5, "Unable to Setup Target");
    ierr = EPSGetST(eps,&st);                                   if (ierr) early_exit(5, "Unable to Setup Target");
    ierr = STSetType(st,"sinvert");                             if (ierr) early_exit(5, "Unable to Setup Target");
    ierr = STSetShift(st, target_eval);                         if (ierr) early_exit(5, "Unable to Setup Target");

    ierr = STGetKSP(st, &ksp);                                  if (ierr) early_exit(6, "Unable to setup KSP");
    ierr = KSPSetType(ksp, KSPPREONLY);                         if (ierr) early_exit(6, "Unable to setup KSP");
    ierr = KSPGetPC(ksp, &pc);                                  if (ierr) early_exit(6, "Unable to setup KSP");
    ierr = PCSetType(pc, PCCHOLESKY);                           if (ierr) early_exit(6, "Unable to setup KSP");

    ierr = EPSSolve(eps);                                       if (ierr) early_exit(7, "Unable to Solve Eigenproblem");

    EPSGetConverged(eps, &nconv);
    if (nconv > 0) {
        ierr = EPSGetEigenpair(eps, 0, &eigenvalue_r, &eigenvalue_i, xr, xi);     if (ierr) early_exit(7, "Unable to get Solution");
        PetscPrintf(PETSC_COMM_WORLD, "eigenvalue: %f", eigenvalue_r);
    } else {
        early_exit(7, "Eigenvalue Solver Failed to converge!");
    }
    EPSDestroy(&eps);

    // ------------------------------------------------------------------------------------------------------------------
    // Return Solution
    // ------------------------------------------------------------------------------------------------------------------
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);            if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = PetscViewerSetType(viewer, "binary");                    if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);         if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = PetscViewerFileSetName(viewer, "./out_vector.dat");      if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = VecView(xr, viewer);                                     if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(8, "Unable to write Eigenvector");


    // ------------------------------------------------------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------------------------------------------------------
    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&xr);
    VecDestroy(&xi);

    SlepcFinalize();
    return 0;

    // // get command line arguments
    // PetscBool flg;
    // PetscScalar target_eval = 1.0;
    // PetscInt dimension = 0;     // square matrix dimensions
    // PetscInt num_values_a = 0;  // number of values in the A Matrix
    // PetscInt num_values_b = 0;  // number of values in the B Matrix
    // char ab_file_names[2][PETSC_MAX_PATH_LEN]; // short unique string associated with this solution (allows multiple non-conflicting executions of this program) 
    
    // PetscOptionsGetScalar(NULL, NULL, "-te", &target_eval, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate target eigenvalue nonce the -te option");
    // PetscOptionsGetInt(NULL, NULL, "-d", &dimension, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate matrix dimension with the -d option");    
    // PetscOptionsGetInt(NULL, NULL, "-va", &num_values_a, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate the number of values in the A matrix with the -va option");   
    // PetscOptionsGetInt(NULL, NULL, "-vb", &num_values_a, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate the number of values in the B matrix with the -vb option");   
    // PetscOptionsGetString(NULL, "-fa", ab_file_names[0], PETSC_MAX_PATH_LEN-1, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate A matrix file name with with the -fa option");
    // PetscOptionsGetString(NULL, "-fb", ab_file_names[1], PETSC_MAX_PATH_LEN-1, &flg);
    // // if (!flg) SETERRQ(2,"Must indicate B matrix file name with with the -fb option");

    // // // load matrix files into SEQAIJ matrices
    // // Mat A, B;
    // // PetscViewer pv_a, pv_b;
    // // PetscViewerBinaryOpen(PETSC_COMM_WORLD,ab_file_names[0],FILE_MODE_READ,&pv_a);
    // // PetscViewerBinaryOpen(PETSC_COMM_WORLD,ab_file_names[1],FILE_MODE_READ,&pv_b);
    // // MatLoad(pv_a,MATSEQAIJ,&A);
    // // MatLoad(pv_b,MATSEQAIJ,&B);
    // // PetscViewerDestroy(pv_a);
    // // PetscViewerDestroy(pv_b);

}