#include "petscviewer.h"
#include <slepceps.h>
#include <petscsys.h>
#include <petscmat.h>
#include <mpi.h>

#include <array>
#include <string>

void early_exit(int status, const char msg[]) {
    PetscPrintf(PETSC_COMM_WORLD, msg);
    SlepcFinalize();
    exit(status);
}

std::array<std::string, 4> file_names(std::string prefix) {
    return {
        prefix + "_a.dat",
        prefix + "_b.dat",
        prefix + "_evec.dat",
        prefix + "_eval.dat",
    };
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
    PetscBool               flg;
    char                    file_prefix[PETSC_MAX_PATH_LEN];

    double eigenvalue_r, eigenvalue_i;
    double target_eval = 1.0;

    // ------------------------------------------------------------------------------------------------------------------
    // Basic Setup and Arguments
    // ------------------------------------------------------------------------------------------------------------------
    int num_threads;
    ierr = SlepcInitialize(&argc, &argv, (char *) 0, NULL); 
    if (ierr) early_exit(1, "Failed to initialize Slepc!");

    MPI_Comm_size(PETSC_COMM_WORLD, &num_threads);
    if (num_threads != 1) early_exit(2, "GEP Solver is only setup to use 1 MPI Thread!");

    PetscOptionsGetScalar(NULL, NULL, "-te", &target_eval, &flg);
    if (!flg) early_exit(2, "Must indicate Target Eigenvalue with the -te option!");

    PetscOptionsGetString(NULL, NULL, "-fp", file_prefix, PETSC_MAX_PATH_LEN-1, &flg);
    if (!flg) early_exit(2, "Must indicate file prefix -fp option!");

    std::array<std::string, 4> filenames = file_names(std::string(file_prefix));

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
    ierr = PetscViewerFileSetName(viewer, &filenames[0][0]);        if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = MatLoad(A, viewer);                                      if (ierr) early_exit(3, "Unable to Create A Matrix");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(3, "Unable to Create A Matrix");

    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);            if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerSetType(viewer, "binary");                    if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);          if (ierr) early_exit(3, "Unable to Create B Matrix");
    ierr = PetscViewerFileSetName(viewer, &filenames[1][0]);        if (ierr) early_exit(3, "Unable to Create B Matrix");
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
    ierr = PetscViewerFileSetName(viewer, &filenames[2][0]);        if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = VecView(xr, viewer);                                     if (ierr) early_exit(8, "Unable to write Eigenvector");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(8, "Unable to write Eigenvector");

    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer);            if (ierr) early_exit(8, "Unable to write Eigenvalue");
    ierr = PetscViewerSetType(viewer, "binary");                    if (ierr) early_exit(8, "Unable to write Eigenvalue");
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);         if (ierr) early_exit(8, "Unable to write Eigenvalue");
    ierr = PetscViewerFileSetName(viewer, &filenames[3][0]);        if (ierr) early_exit(8, "Unable to write Eigenvalue");
    ierr = PetscScalarView(1, &eigenvalue_r, viewer);               if (ierr) early_exit(8, "Unable to write Eigenvalue");
    ierr = PetscViewerDestroy(&viewer);                             if (ierr) early_exit(8, "Unable to write Eigenvalue");

    // ------------------------------------------------------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------------------------------------------------------
    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&xr);
    VecDestroy(&xi);

    SlepcFinalize();
    return 0;
}