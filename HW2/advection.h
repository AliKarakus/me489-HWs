
#ifndef ADVECTION_H
#define ADVECTION_H 1
  # include <stdio.h>
  # include <stdlib.h>
  # include <math.h>
  # include <time.h>
  # include<string.h>  

typedef struct{
    int Nnodes, NX, NY; 
    double xmin, xmax, ymin, ymax; 
    double *x, *y; 
    // Node connectivity 
    int *N2N; 
  }mesh_t;

  typedef struct{
    int Nstage; 
    double time, tstart, tend, dt; 
    double *rhsq, *resq; 
    double *rk4a, *rk4b, *rk4c; 
  }tstep_t;

typedef struct{
    mesh_t msh; 
    tstep_t tstep; 
    double *q, *u; 
  }solver_t;

  // Function definitions
  mesh_t createMesh(char *inputFile); 
  tstep_t createTimeStepper(int Nnodes);
  void initialCondition(solver_t *solver);
  void RhsQ(solver_t *solver, tstep_t *tstep, int stage);
  void solverPlot(char *fileName, mesh_t *msh, double *Q);
  double readInputFile(char *fileName, char *tag);
 

#endif