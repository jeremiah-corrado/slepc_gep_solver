# slepc_gep_solver

A Sparse Generalized Eigenvalue Solver intended to be used with the [fem_2d library](https://github.com/jeremiah-corrado/fem_2d)

## Depedancies

* [Make](https://cmake.org/)
* [MPICH](https://www.mpich.org/)
* [Petsc](https://petsc.org/release/overview/)
* [Slepc](https://slepc.upv.es/)

Linux is highly recomended, but it should be possible to install Slepc and Petsc on Windows or Mac if needed.

## Instalation

1. Install Petsc. You can follow their instructions [here](https://petsc.org/release/install/install_tutorial/)
2. Install Slepc (an extension of Petsc). You can follow their instructions in the installation section of [this](https://petsc.org/release/install/install_tutorial/) pdf
3. Download this repository into your directory of choice
4. Navigate to your dicrectory in a terminal and run `make`. This should produce a binary file called **solve_gep**
5. Set the environment variable `GEP_SOLVE_DIR` to the directory where **solve_gep** and the **/tmp/** folder are located. This will allow *fem_2d* to find the solver.
  * This variable can be set perminantly by adding the line GEP_SOLVE_DIR="the/actual/path/on/your/system/slepc_gep_solver/" to your `.environment` file in `/etc/`
