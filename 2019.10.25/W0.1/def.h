#include <petscdmda.h>
#include <petsctime.h>
#include <math.h>

#define EXAMPLE 2               /* 2 for Poiseuille flow */
#define CHECKERROR 1            /* 1 for question that have Getanalysis */
#define highorderext 1          /* 1 for second order accuracy NEM */
#define PI 3.14159265354

#if EXAMPLE == 2                /* Poiseuille flow */
    #define BOUNDARYX 1         /* periodic */
    #define BOUNDARYY 0         /* noperiodic */
#endif

typedef struct {
  PetscScalar    f[9];
} Field;

typedef struct {
  PetscReal      lidvelocity, Re, nu, rho0;             /* physical parameters */
  PetscBool      draw_contours;                /* flag - 1 indicates drawing contours */
  PetscScalar    eu[9], ev[9];
  PetscScalar    w[9], W, rate;
  PetscInt       c, nits;
  PetscScalar    c0;                           /* dx/dt */
  PetscScalar    crtime, ptime, stime, ee[2000][3], step, F[200],tao[200];              /* current, previous time and time step Modified by program */
  Vec            x, px, func;                              /* current and previous solution  */
  PetscInt       Nx, Q, Ny;
  PetscScalar    *error, Linfty, L1, L2, maxtime, ep;
} AppCtx;

extern PetscErrorCode InitialGuess(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode evolution(DMDALocalInfo*, Vec, AppCtx*);
extern PetscScalar Error(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode SetOptions(AppCtx*);
extern PetscScalar CalFeq(PetscInt, PetscScalar, PetscScalar, PetscScalar, AppCtx*);
extern PetscErrorCode DataSaveASCII(DMDALocalInfo*, Vec, char*, char*);
extern PetscErrorCode DataSaveASCIIF(DMDALocalInfo*, AppCtx*, char*);
extern PetscErrorCode Getanalysis(DMDALocalInfo*, AppCtx*);
extern PetscErrorCode ComputeLfunction(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode ComputeLBoundary(DMDALocalInfo*, Field**, AppCtx*);
extern PetscErrorCode InitialF(DMDALocalInfo*, AppCtx*);
extern PetscErrorCode InitialTAO(DMDALocalInfo*, AppCtx*);