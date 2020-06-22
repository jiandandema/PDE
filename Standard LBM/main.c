
static char help[] = "/********* The Poiseuille flow in 2d *************/\n\
  \n\
The Poiseuille flow is a channel flow driven by a constant force along the x1\n\
direction between two parallel plates. There is a steady state solution that\n\
can be expressed as a parabola centered around the axis of the channel.\n\
this numerical experiments refer to section 4.1 of paper 'A FULLY IMPLICIT\n\
METHOD FOR LATTICE BOLTZMANNEQUATIONS'\n\
  -lidvelocity <lid>, where <lid> = dimensionless velocity of lid\n\
  -contours : draw contour plots of solution\n"; 

#include <petscdmda.h>
#include "def.h"

int main(int argc, char **argv)
{
  AppCtx         user; 
  PetscInt       mx, my, xs, ys, xm, ym, N, MaxN, nits; 
  DMDALocalInfo  info; 
  MPI_Comm       comm; 
  DM             da;
  DMBoundaryType bx = DM_BOUNDARY_PERIODIC, by = DM_BOUNDARY_NONE; 
  PetscScalar    error = 1; 
  char           sn[40], sm[40];
  PetscLogDouble v2, vv1;                /* Returns the CPU time in seconds used by the process */

  PetscInitialize(&argc, &argv, (char*)0, help); 
  comm = PETSC_COMM_WORLD; 
  N = 101; 
  MaxN = 1000000; 

  PetscPreLoadBegin(PETSC_FALSE, "SetUp");

  DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR, N, N, PETSC_DECIDE, PETSC_DECIDE, 9, 2, NULL, NULL, &da); 
  DMSetFromOptions(da); 
  DMSetUp(da); 
  DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);  
  DMSetApplicationContext(da, &user); 
  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, 
              PETSC_IGNORE, PETSC_IGNORE); 
  DMDAGetLocalInfo(da, &info); 
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); 
    
  user.Nx = mx; 
  SetOptions(&user); 
  PetscOptionsGetReal(NULL, NULL, "-lidvelocity", &user.lidvelocity, NULL); 
  PetscOptionsGetReal(NULL, NULL, "-Re", &user.Re, NULL);                /* Gets the double precision value for a particular option in the database */

  DMDASetFieldName(da, 0, "f0"); 
  DMDASetFieldName(da, 1, "f1"); 
  DMDASetFieldName(da, 2, "f2"); 
  DMDASetFieldName(da, 3, "f3"); 
  DMDASetFieldName(da, 4, "f4"); 
  DMDASetFieldName(da, 5, "f5"); 
  DMDASetFieldName(da, 6, "f6"); 
  DMDASetFieldName(da, 7, "f7"); 
  DMDASetFieldName(da, 8, "f8"); 

  DMCreateGlobalVector(da, &user.x); 
  DMCreateGlobalVector(da, &user.func); 
  DMDestroy(&da); 
  
  InitialGuess(&info, user.x, &user); 
  PetscTime(&vv1);
  for (nits = 1; nits <= MaxN; nits++) {
    user.ptime = user.crtime; 
    user.crtime += user.stime; 
    user.nits = nits; 
#if CHECKERROR == 1
    Getanalysis(&info, &user); 
#endif	
    if (nits == 1) {
      VecDuplicate(user.x, &user.px); 
      VecCopy(user.x, user.px); 
    }
    if (nits > 1) {
      VecDestroy(&user.px); 
      VecDuplicate(user.x, &user.px); 
      VecCopy(user.x, user.px); 
    }

    evolution(&info, user.x, &user); 


    if (nits % 1000 == 1) {
      PetscPrintf(comm, "\n/********* start new time layer %d *************/\n", nits); 
      PetscPrintf(comm, "current time layer: [%g, %g]\n", user.ptime, user.crtime); 
      sprintf(sn, "pre_ex%dout.dat", nits);                /* Prints information to a PetscViewer string */
      sprintf(sm, "ex%dout.dat", nits);   
      DataSaveASCII(&info, user.x, sn, sm); 
#if CHECKERROR == 1
      error = Error(&info, user.x, &user); 
      PetscPrintf(comm, "current error: [%g]\n", error);
#endif
    }

    if (error < 1e-6 || nits == MaxN || user.crtime + 1e-6 > user.maxtime) {
      PetscPrintf(comm, "current time layer: [%g, %g]\n", user.ptime, user.crtime);       
      sprintf(sn, "pre_exout.dat");                /* Prints information to a PetscViewer string */
      sprintf(sm, "exout.dat");   
      DataSaveASCII(&info, user.x, sn, sm); 
      break; 
    }
  }
  PetscTime(&v2);
  
  PetscPrintf(PETSC_COMM_WORLD, "Total CPU time = %g\n", v2 - vv1); 
  PetscPrintf(PETSC_COMM_WORLD, "Last compute time = %g\n", user.L1); 

  PetscPreLoadEnd(); 
  PetscFinalize();   
  return 0; 
}


PetscErrorCode InitialGuess(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt       i, j, k; 
  Field          **x; 
  PetscScalar	 rho, u, v; 

  DMDAVecGetArray(info->da, X, &x); 

  //p0 = 1.0;                /* the initialguess of th0 and p0 must be coupled with Poisson equation */
  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      rho = 1.0; u = 0.0; v = 0.0; 
      for (k = 0; k < user->Q; k++) {
        x[j][i].f[k]  = CalFeq(k, rho, u, v, user); 
      }
    }
  }

  DMDAVecRestoreArray(info->da, X, &x); 
  return 0; 
}


PetscErrorCode evolution(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  Field 	 **x, **xp; 
  Vec		 loc_xp;
  PetscLogDouble v1, v2, v3; 

  PetscFunctionBegin;                /* First executable line of each PETSc function used for error handling */

  DMGetLocalVector(info->da, &loc_xp); 
  DMGlobalToLocalBegin(info->da, user->px, INSERT_VALUES, loc_xp);
  DMGlobalToLocalEnd(info->da, user->px, INSERT_VALUES, loc_xp);
  DMDAVecGetArray(info->da, loc_xp, &xp);                /* the local xp[j][i] aboout local px = x */
  DMDAVecGetArray(info->da, X, &x);                /* the global x[j][i] aboout global x */

  PetscTime(&v1);


  Collision_Step(info, xp, x, user);
  Propagation_Step(info, x, user);
  ComputeLBoundary(info, x, user);

  PetscTime(&v2);
  v3 = v2 - v1;
  user->L1 = v3;

  DMDAVecRestoreArray(info->da, X, &x); 
  DMDAVecRestoreArray(info->da, loc_xp, &xp); 
  VecDestroy(&loc_xp); 

  PetscFunctionReturn(0); 
}


PetscErrorCode SetOptions(AppCtx *user)
{
  PetscInt       c = 1; 
    
  user->c = 1.0; user->c0 = 1.0; 
  user->eu[0] = c*0; 	  user->ev[0] = 0*c;    user->w[0]=4.0/9.0; 
  user->eu[1] = c*1;    user->ev[1] = 0*c;    user->w[1]=1.0/9.0; 
  user->eu[2] = c*0;    user->ev[2] = 1*c;    user->w[2]=1.0/9.0; 
  user->eu[3] = -1*c;   user->ev[3] = 0*c;    user->w[3]=1.0/9.0; 
  user->eu[4] = 0*c;    user->ev[4] = -1*c;   user->w[4]=1.0/9.0; 
  user->eu[5] = 1*c;    user->ev[5] = 1*c;    user->w[5]=1.0/36.0; 
  user->eu[6] = -1*c;   user->ev[6] = 1*c;    user->w[6]=1.0/36.0; 
  user->eu[7] = -1*c;   user->ev[7] = -1*c;   user->w[7]=1.0/36.0; 
  user->eu[8] = 1*c;    user->ev[8] = -1*c;   user->w[8]=1.0/36.0; 

  user->crtime = 0.0;
  user->ptime = 0.0;
  user->H_star = 100.0;
  user->rho_star = 1.0;
  user->tau_star = 0.55;
  //user->nu = 1e-6;
  //user->step = 1e-5;
  user->stime =  1e-4 * 5.0 / 3.0;                /* for CFL number be small */
  user->maxtime = 200;
  user->lidvelocity = 1.25 / 6.0;
  user->Q = 9;

  return 0; 
}


PetscScalar CalFeq(PetscInt k, PetscScalar rho, PetscScalar u, PetscScalar v, AppCtx *user)
{
  PetscScalar    eu, uv, feq, c; 
  
  c = user->c0; 
  eu = (c * user->eu[k] * u + c * user->ev[k] * v) / (PetscScalar)(user->c) / (PetscScalar)(user->c); 
  uv = (u * u + v * v ) / (PetscScalar)(user->c) / (PetscScalar)(user->c); 
  feq = user->w[k] * rho * (1 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv); 
  
  return feq; 
}
