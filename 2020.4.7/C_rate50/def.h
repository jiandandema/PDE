#include <petscdmda.h>
#include <petsctime.h>
#include <math.h>
#define pi 3.14159265358979323846

#define EXAMPLE 2      /* 2 for Poiseuille flow */
#define CHECKERROR 1   /* 1 for question that have Getanalysis */
#define highorderext 1 /* 1 for second order accuracy NEM */

#if EXAMPLE == 2    /* Poiseuille flow */
#define BOUNDARYX 1 /* periodic */
#define BOUNDARYY 0 /* noperiodic */
#endif

typedef struct
{
  PetscScalar f[9];
} Field;

typedef struct
{
  PetscReal lidvelocity, rho0; /* physical parameters */
  PetscBool draw_contours;     /* flag - 1 indicates drawing contours */
  PetscScalar eu[9], ev[9], W, rate, H, nu_l;
  PetscScalar w[9];
  PetscInt c, nits;
  PetscScalar c0;                                      /* dx/dt */
  PetscScalar crtime, ptime, stime, ee[2000][3], step; /* current, previous time and time step Modified by program */
  Vec x, px, func, tao;                                /* current and previous solution  */
  PetscInt Nx, Q, Ny;
  PetscScalar *error, Linfty, L1, L2, maxtime, ep;
} AppCtx;

extern PetscErrorCode InitialGuess(DMDALocalInfo *, Vec, AppCtx *);
extern PetscErrorCode InitialTao(DMDALocalInfo *, Vec, AppCtx *);
extern PetscErrorCode InitialF(DMDALocalInfo *, Vec, AppCtx *);
extern PetscErrorCode evolution(DMDALocalInfo *, Vec, Vec, AppCtx *);
extern PetscScalar Error(DMDALocalInfo *, Vec, AppCtx *);
extern PetscErrorCode SetOptions(AppCtx *);
extern PetscScalar CalFeq(PetscInt, PetscScalar, PetscScalar, PetscScalar, AppCtx *);
extern PetscErrorCode DataSaveASCII(DMDALocalInfo *, Vec, char *, char *);
extern PetscErrorCode Getanalysis(DMDALocalInfo *, AppCtx *);
extern PetscErrorCode ComputeLfunction(DMDALocalInfo *, Field **, Field **, Field **, AppCtx *);
extern PetscErrorCode ComputeLBoundary(DMDALocalInfo *, Field **, AppCtx *);
extern PetscErrorCode DataSaveASCIITao(DMDALocalInfo *, Vec, char *, char *);