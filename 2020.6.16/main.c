
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
  PetscScalar    error = 1, perror=1, eerror = 1; 
  char           sn[40], sm[40];
  PetscLogDouble v2, vv1;                /* Returns the CPU time in seconds used by the process */

  PetscInitialize(&argc, &argv, (char*)0, help); 
  comm = PETSC_COMM_WORLD; 
  N = 17; 
  MaxN = 1000000; 

  PetscPreLoadBegin(PETSC_FALSE, "SetUp");

  DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, N, N, PETSC_DECIDE, PETSC_DECIDE, 9, 2, NULL, NULL, &da); 
  DMSetFromOptions(da); 
  DMSetUp(da); 
  DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);  
  DMSetApplicationContext(da, &user); 
  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, 
              PETSC_IGNORE, PETSC_IGNORE); 
  DMDAGetLocalInfo(da, &info); 
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); 
    
  user.Nx = mx; 
  user.Ny = my;
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
  DMCreateGlobalVector(da, &user.f_neq); 
  DMCreateGlobalVector(da, &user.nu);
  DMCreateGlobalVector(da, &user.grad_nu);
  DMCreateGlobalVector(da, &user.F);
  DMDestroy(&da); 
  
  InitialGuess(&info, user.x, &user);

  InitialNu(&info, user.nu, &user);
  sprintf(sn, "pre_exnu.dat");                /* Prints information to a PetscViewer string */
  sprintf(sm, "exnu.dat");
  DataSaveASCII(&info, user.nu, sn, sm);

  InitialGradNu(&info, user.grad_nu, &user);
  sprintf(sn, "pre_exgrad_nu.dat");                /* Prints information to a PetscViewer string */
  sprintf(sm, "exgrad_nu.dat");
  DataSaveASCII(&info, user.grad_nu, sn, sm);
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
    InitialFneq(&info, user.x, user.f_neq, &user);
    InitialF(&info, user.F, user.f_neq, user.nu, user.grad_nu, &user);
    evolution(&info, user.x, user.F, &user); 

    if (nits % 1000 == 1) {
      PetscPrintf(comm, "\n/********* start new time layer %d *************/\n", nits); 
      PetscPrintf(comm, "current time layer: [%g, %g]\n", user.ptime, user.crtime); 
      sprintf(sn, "pre_ex%dout.dat", nits);                /* Prints information to a PetscViewer string */
      sprintf(sm, "ex%dout.dat", nits);   
      DataSaveASCII(&info, user.x, sn, sm); 
#if CHECKERROR == 1
      error = Error(&info, user.x, &user);
      eerror = fabs(error - perror);
      perror = error; 
      PetscPrintf(comm, "current error: [%g]\n", error);
#endif
    }

    if (eerror < 1e-6 || nits == MaxN || user.crtime + 1e-6 > user.maxtime) {
      PetscPrintf(comm, "current time layer: [%g, %g]\n", user.ptime, user.crtime);       
      sprintf(sn, "pre_exout.dat");                /* Prints information to a PetscViewer string */
      sprintf(sm, "exout.dat");   
      DataSaveASCII(&info, user.x, sn, sm);
      sprintf(sn, "pre_exF.dat");                /* Prints information to a PetscViewer string */
      sprintf(sm, "exF.dat");
      DataSaveASCII(&info, user.F, sn, sm);
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


PetscErrorCode InitialFneq(DMDALocalInfo *info, Vec X, Vec Y, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k;
  Field          **f_neq, **x;
  PetscScalar    sumrho, sumu, sumv;

  PetscFunctionBegin;
  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym;

  DMDAVecGetArray(info->da, X, &x);
  DMDAVecGetArray(info->da, Y, &f_neq);

  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      sumrho = 0.0; sumu = 0.0; sumv = 0.0;
      for (k = 0; k < user->Q; k++) {
          sumrho += x[j][i].f[k];
          sumu   += user->eu[k] * x[j][i].f[k];
          sumv   += user->ev[k] * x[j][i].f[k];
        }
        sumu = sumu / sumrho;
        sumv = sumv / sumrho;
      for (k = 0; k < user->Q; k++)
        f_neq[j][i].f[k] = x[j][i].f[k] - CalFeq(k, sumrho ,sumu, sumv, user);
    }
  }

  DMDAVecRestoreArray(info->da, Y, &f_neq);
  DMDAVecRestoreArray(info->da, X, &x);
  PetscFunctionReturn(0);
}



PetscErrorCode InitialNu(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt i, j;
  Field **x;
  PetscScalar nu = 0;
  PetscScalar phi = 0, n = 0;

  DMDAVecGetArray(info->da, X, &x);

  for (j = info->ys; j < info->ys + info->ym; j++)
  {
    for (i = info->xs; i < info->xs + info->xm; i++)
    {
      if (j < user->Ny / 2.0)
      {
        n = tanh((j * user->step - 1.0 / 3.0) / user->W);
        phi = 0.5 * (-n + 1);
        nu = ((user->rate - 1) * phi + 1) * user->eta_l;
        x[j][i].f[0] = nu;
      }
      else
      {
        n = tanh((j * user->step - 2.0 / 3.0) / user->W);
        phi = 0.5 * (n + 1);
        nu = ((user->rate - 1) * phi + 1) * user->eta_l;
        x[j][i].f[0] = nu;
      }
    }
  }

  DMDAVecRestoreArray(info->da, X, &x);
  return 0;
}



PetscErrorCode InitialGradNu(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt i, j;
  PetscScalar n = 0;
  Field **x;

  DMDAVecGetArray(info->da, X, &x);

  for (j = info->ys; j < info->ys + info->ym; j++)
  {
    for (i = info->xs; i < info->xs + info->xm; i++)
    {
      if (j < user->Ny / 2.0)
      {
       n = tanh((j * user->step - 1.0 / 3.0) / user->W);
       x[j][i].f[1] = -0.5 * (user->rate - 1) * user->eta_l * (1 - n * n) / user->W;
       x[j][i].f[0] = 0;
      }
      else
      {
        n = tanh((j * user->step - 2.0 / 3.0) / user->W);
        x[j][i].f[1] = 0.5 * (user->rate - 1) * user->eta_l * (1 - n * n) / user->W;
        x[j][i].f[0] = 0;
      }
    }
  }

  DMDAVecRestoreArray(info->da, X, &x);
  return 0;
}


PetscErrorCode InitialF(DMDALocalInfo *info, Vec X, Vec Y, Vec Z, Vec W, AppCtx *user)
{
   PetscInt       xints, xinte, yints, yinte, i, j, k, l;
   Field          **x, **f_neq, **nu, **grad_nu, **u;
   PetscScalar    S_xx = 0, S_yy = 0, S_xy = 0, S_yx = 0, S_xx_x = 0, S_yy_y = 0, S_xy_y = 0, S_yx_x = 0;
   float          S_xx_dev[9], S_xy_dev[9], S_yy_dev[9], S_yx_dev[9];
   Vec            loc_f_neq, loc_x;

  PetscFunctionBegin;
  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym;

  DMDAVecGetArray(info->da, X, &x);
  DMGetLocalVector(info->da, &loc_f_neq); 
  DMGlobalToLocalBegin(info->da, Y, INSERT_VALUES, loc_f_neq); 
  DMGlobalToLocalEnd(info->da, Y, INSERT_VALUES, loc_f_neq); 
  DMDAVecGetArray(info->da, loc_f_neq, &f_neq);
  DMDAVecGetArray(info->da, Z, &nu);
  DMDAVecGetArray(info->da, W, &grad_nu);


  DMGetLocalVector(info->da, &loc_x); 
  DMGlobalToLocalBegin(info->da, user->x, INSERT_VALUES, loc_x); 
  DMGlobalToLocalEnd(info->da, user->x, INSERT_VALUES, loc_x); 
  DMDAVecGetArray(info->da, loc_x, &u);


  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {


        S_xx = 0;S_yy = 0;S_xy = 0;S_yx = 0;S_xx_x = 0;S_yx_x = 0;S_xy_y = 0;S_yy_y = 0;
        for (k = 0; k < user->Q; k++)
        {
          S_xx_dev[k] = 0;S_xy_dev[k] = 0;S_yx_dev[k] = 0;S_yy_dev[k] = 0;
          x[j][i].f[k] = 0;
        }


        for(k = 0; k < user->Q; k++)
        {
          S_xx += user->eu[k] * user->eu[k] * f_neq[j][i].f[k];
          S_yy += user->ev[k] * user->ev[k] * f_neq[j][i].f[k];
          S_xy += user->eu[k] * user->ev[k] * f_neq[j][i].f[k];
          S_yx += user->ev[k] * user->eu[k] * f_neq[j][i].f[k];
          for(l = 0; l < user->Q; l++)
          {
            S_xx_dev[k] += user->eu[l] * user->eu[l] * f_neq[j + user->ev[k]][i + user->eu[k]].f[l];
            S_xy_dev[k] += user->eu[l] * user->ev[l] * f_neq[j + user->ev[k]][i + user->eu[k]].f[l];
            S_yx_dev[k] += user->ev[l] * user->eu[l] * f_neq[j + user->ev[k]][i + user->eu[k]].f[l];
            S_yy_dev[k] += user->ev[l] * user->ev[l] * f_neq[j + user->ev[k]][i + user->eu[k]].f[l];
          }

          S_xx_x += user->eu[k] * S_xx_dev[k];
          S_xy_y += user->ev[k] * S_xy_dev[k];
          S_yx_x += user->eu[k] * S_yx_dev[k];
          S_yy_y += user->ev[k] * S_yy_dev[k];
        }
        S_xx_x = (user->Nx - 1) / 6.0 * S_xx_x;
        S_yx_x = (user->Nx - 1) / 6.0 * S_yx_x;
        S_xy_y = (user->Nx - 1) / 6.0 * S_xy_y;
        S_yy_y = (user->Nx - 1) / 6.0 * S_yy_y;
        S_xx = -0.5 / user->eta_l * S_xx;
        S_xy = -0.5 / user->eta_l * S_xy;
        S_yx = -0.5 / user->eta_l * S_yx;
        S_yy = -0.5 / user->eta_l * S_yy;
        S_xx_x = -0.5 / user->eta_l * S_xx_x;
        S_yx_x = -0.5 / user->eta_l * S_yx_x;
        S_xy_y = -0.5 / user->eta_l * S_xy_y;
        S_yy_y = -0.5 / user->eta_l * S_yy_y;


        //x[j][i].f[0] = 2.0 * S_xx * grad_nu[j][i].f[0] + 2.0 * S_xy * grad_nu[j][i].f[1] + 2.0 * (nu[j][i].f[0] - user->eta_l) * S_xx_x + 2.0 * (nu[j][i].f[0] - user->eta_l) * S_xy_y;
        x[j][i].f[0] = 2.0 * S_xy * grad_nu[j][i].f[1] + 2.0 * (nu[j][i].f[0] - user->eta_l) * S_xy_y;
        //x[j][i].f[1] = 2.0 * S_yy * grad_nu[j][i].f[1] + 2.0 * S_yx * grad_nu[j][i].f[0] + 2.0 * (nu[j][i].f[0] - user->eta_l) * S_yx_x + 2.0 * (nu[j][i].f[0] - user->eta_l) * S_yy_y;
        x[j][i].f[1] = 0;
        x[j][i].f[0] += 8 * user->lidvelocity * user->rho0 * user->rate * user->eta_l;
        //x[j][i].f[0] += 0;
      }
    }
  }


  DMDAVecRestoreArray(info->da, loc_f_neq, &f_neq);
  DMDAVecRestoreArray(info->da, Z, &nu);
  DMDAVecRestoreArray(info->da, W, &grad_nu);
  DMDAVecRestoreArray(info->da, X, &x);
  DMDAVecRestoreArray(info->da, loc_x, &u);
  VecDestroy(&loc_f_neq);
  VecDestroy(&loc_x);
  PetscFunctionReturn(0);
}


PetscErrorCode evolution(DMDALocalInfo *info, Vec X, Vec F, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  Field 	 **x, **xp, **x1, **loc_X1, **f; 
  Vec		 loc_x1, loc_xp, X1;
  PetscLogDouble v1, v2, v3; 

  PetscFunctionBegin;                /* First executable line of each PETSc function used for error handling */

  DMGetLocalVector(info->da, &loc_xp); 
  DMGlobalToLocalBegin(info->da, user->px, INSERT_VALUES, loc_xp); 
  DMGlobalToLocalEnd(info->da, user->px, INSERT_VALUES, loc_xp); 
  DMDAVecGetArray(info->da, loc_xp, &xp);                /* the local xp[j][i] aboout local px = x */
  VecDuplicate(X, &X1); 
  DMDAVecGetArray(info->da, X1, &x1);                /* the global x1[j][i] with no value global X1 */
  DMDAVecGetArray(info->da, X, &x);                /* the global x[j][i] aboout global x */
  DMDAVecGetArray(info->da, F, &f);

  PetscTime(&v1);
  ComputeLfunction(info, xp, x1, f, user);                /* compute inner point for x1 */
  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym;                /* same info from DMDALocalInfo *info */
  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {
        for (k = 0; k < user->Q; k++) {
          x1[j][i].f[k] = xp[j][i].f[k] + user->stime * x1[j][i].f[k];                /* refer to ComputeLfunction() */
        }
      }
    }
  }
  ComputeLBoundary(info, x1, user);                /* compute boundry point for x1 */
  
  PetscTime(&v2);
  v3 = v2 - v1;
  DMDAVecRestoreArray(info->da, X1, &x1);                /* X1 is unrelated to x1[j][i]  */

  DMGetLocalVector(info->da, &loc_x1); 
  DMGlobalToLocalBegin(info->da, X1, INSERT_VALUES, loc_x1); 
  DMGlobalToLocalEnd(info->da, X1, INSERT_VALUES, loc_x1); 
  DMDAVecGetArray(info->da, loc_x1, &loc_X1);                /* the local loc_X1[j][i] aboout local X1 */

  PetscTime(&v1);
  ComputeLfunction(info, loc_X1, x, f, user); 
  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {
        for (k = 0; k < user->Q; k++) {
          x[j][i].f[k] = 0.5 * (xp[j][i].f[k] + loc_X1[j][i].f[k]) + 0.5 * user->stime * x[j][i].f[k];                /* second order accuracy Runge-Kutta scheme (SSP RK-2)*/
        }
      }
    }
  }
  ComputeLBoundary(info, x, user);
  
  PetscTime(&v2);
  v3 += v2 - v1;
  user->L1 = v3;

  DMDAVecRestoreArray(info->da, X, &x); 
  DMDAVecRestoreArray(info->da, F, &f);
  DMDAVecRestoreArray(info->da, loc_xp, &xp); 
  DMDAVecRestoreArray(info->da, loc_x1, &loc_X1); 
  VecDestroy(&X1); 
  VecDestroy(&loc_xp); 
  VecDestroy(&loc_x1); 
  
  PetscFunctionReturn(0); 
}


PetscErrorCode SetOptions(AppCtx *user)
{
  PetscInt       c = 1; 
    
  user->c = 1; user->c0 = 1.0; 
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
  user->step = 1.0 / (user->Nx - 1);
  user->ptime = 0.0; 
  user->stime = 0.1 / (user->Nx - 1);                /* for CFL number be small */
  user->rho0 = 1.0; 
  user->Re  = 0.5; 
  user->maxtime = 200;
  user->lidvelocity = 0.001;
  user->eta_l = 0.002;                /* viscosity coefficent */
  user->tau = 3.0 * user->eta_l;                /* maybe the one oder accuracy scheme for dimensionless tau and viscosity coefficent */
  user->W = 0.01;
  user->rate = 50;
  //user->tau = 3.0 * user->nu * 0.01;                /* 使得碰撞项占优而不是F占优，同时也会放大迎风格式带来的非物理的扰动（比如ep=0.1或0.0时） */
  //user->tau = 3.0 * user->nu + 0.5 * user->stime;                /* I think &t in LBM should be user->stime */
  //user->tau = 3.0 * user->nu + 0.5 / (user->Nx - 1);                 /* I think &t in LBM may be 1 / (user->Nx - 1) */
  //user->F_x = 8 * user->lidvelocity * user->rho0 * user->nu;
  //user->F_y = 0;
  user->Q = 9; 
  user->ep = 0.0; 
  PetscOptionsGetReal(NULL, NULL, "-ep", &user->ep, NULL); 

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
