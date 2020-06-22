#include <petscdmda.h>
#include <petsctime.h>

#define EXAMPLE 2               /* 2 for Poiseuille flow */
#define CHECKERROR 1            /* 1 for question that have Getanalysis */ 
#define highorderext 1          /* 1 for second order accuracy NEM */  

#if EXAMPLE == 2                /* Poiseuille flow */
    #define BOUNDARYX 1         /* periodic */
    #define BOUNDARYY 0         /* noperiodic */
#endif

typedef struct {
  PetscScalar    f[9];
} Field;

typedef struct {
  PetscReal      lidvelocity, Re, rho0, tau;             /* physical parameters */
  PetscBool      draw_contours;                /* flag - 1 indicates drawing contours */
  PetscScalar    rate, eta_l, W;
  PetscScalar    w[9];
  PetscInt       c, nits, eu[9], ev[9];
  PetscScalar    c0, step;                           /* dx/dt */
  PetscScalar    crtime, ptime, stime, ee[2000][3];              /* current, previous time and time step Modified by program */
  Vec            x, px, func, S, nu, grad_nu, grad_s;                              /* current and previous solution  */
  PetscInt       Nx, Ny, Q;
  PetscScalar    *error, Linfty, L1, L2, maxtime, ep;
} AppCtx;

extern PetscErrorCode InitialGuess(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode evolution(DMDALocalInfo*, Vec, AppCtx*);
extern PetscScalar Error(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode SetOptions(AppCtx*);
extern PetscScalar CalFeq(PetscInt, PetscScalar, PetscScalar, PetscScalar, AppCtx*);
extern PetscErrorCode DataSaveASCII(DMDALocalInfo*, Vec, char*, char*);
extern PetscErrorCode Getanalysis(DMDALocalInfo*, AppCtx*);
extern PetscErrorCode ComputeLfunction(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode ComputeLBoundary(DMDALocalInfo*, Field**, AppCtx*);
//extern PetscErrorCode InitialF(DMDALocalInfo*, Vec, Vec, Vec, Vec, AppCtx*);
extern PetscErrorCode InitialGradNu(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode InitialNu(DMDALocalInfo*, Vec, AppCtx*);
extern PetscErrorCode InitialS(DMDALocalInfo*, Vec, Vec, AppCtx*);
extern PetscErrorCode InitialGradS(DMDALocalInfo*, Vec, Vec, AppCtx*);