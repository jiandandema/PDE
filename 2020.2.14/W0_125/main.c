
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
  N = 200; 
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
  user.Ny = my;
  user.step = 1.0 / (user.Ny - 1);
  SetOptions(&user); 
  PetscOptionsGetReal(NULL, NULL, "-lidvelocity", &user.lidvelocity, NULL); 
  //PetscOptionsGetReal(NULL, NULL, "-Re", &user.Re, NULL);                /* Gets the double precision value for a particular option in the database */

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
  DMCreateGlobalVector(da, &user.tao);
  DMCreateGlobalVector(da, &user.func);
  DMDestroy(&da); 
  
  InitialTao(&info, user.tao, &user);
  sprintf(sn, "pre_tao.dat");
  sprintf(sm, "extao.dat");
  DataSaveASCIITao(&info, user.tao, sn, sm); 

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

    evolution(&info, user.x, user.tao, &user); 

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


PetscErrorCode InitialTao(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt       i, j, k; 
  Field          **x; 
  PetscScalar	 nu; 
  PetscScalar  phi, n;

  DMDAVecGetArray(info->da, X, &x); 

  for (j = info->ys; j < info->ys + info->ym; j++) {
        for (i = info->xs; i < info->xs + info->xm; i++) {
                if (j < user->Ny / 2) {
                      n = tanh((j * user->step - 1.0 / 3.0) / user->W);
                        phi = 0.5 * (-n + 1);
                        nu =  user->rate * phi + 1;
                        x[j][i].f[0] = nu * 3.0 / user->Re;
		}
                else {
                        n = tanh((j * user->step - 2.0 / 3.0) / user->W);
                        phi = 0.5 * (n + 1);
                        nu = user->rate * phi + 1;
                        x[j][i].f[0] = nu * 3.0 / user->Re;
                }
                for (k = 1; k < user->Q; k++) {
                        x[j][i].f[k]  = x[j][i].f[0];
                }
        }
  }

  DMDAVecRestoreArray(info->da, X, &x);
  return 0; 
}


/*
PetscErrorCode InitialF(DMDALocalInfo *info, Vec X, AppCtx *user){
  PetscInt        i,j,k;
  Field           **x;
  PetscScalar     grad_nu, grad_u, grad_r;
  PetscScalar     n, r, xl, yl;

  DMDAVecGetArray(info->da, X, &x);
  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      xl = i*user->step;yl = j*user->step;
      r = sqrt(pow(xl - 1.0/2.0,2)+pow(yl - 1.0/2.0,2));
      n = tanh(2.4*(r - user->H / 16.0) / user->W);
      grad_r = (j * user->step - 1.0/2.0)/r;
      n = 1 - n*n;
      grad_nu =  -(user->rate - 1.0)/(2.0 * user->W) * n * grad_r * 2.4;
      grad_u = user->lidvelocity;
      x[j][i].f[0] = -grad_nu * grad_u;
      x[j][i].f[0] = x[j][i].f[0] / user->Re;
      for (k = 1; k < user->Q; k++) {
        x[j][i].f[k]  = x[j][i].f[0];
      }
    }
  }
  DMDAVecRestoreArray(info->da, X, &x);
  return 0;
}
*/


PetscErrorCode InitialGuess(DMDALocalInfo* info, Vec X, AppCtx* user)
{
	PetscInt       i, j, k;
	Field 		 **x;
	PetscScalar	 rho, u, v;

	DMDAVecGetArray(info->da, X, &x);

	//p0 = 1.0;                /* the initialguess of th0 and p0 must be coupled with Poisson equation */
	for (j = info->ys; j < info->ys + info->ym; j++) {
		for (i = info->xs; i < info->xs + info->xm; i++) {
			rho = 1.0; u = 0.0; v = 0.0;
			for (k = 0; k < user->Q; k++) {
				x[j][i].f[k] = CalFeq(k, rho, u, v, user);
			}
		}
	}

	DMDAVecRestoreArray(info->da, X, &x);
	return 0;
}


PetscErrorCode evolution(DMDALocalInfo *info, Vec X, Vec Y, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  Field 	 **x, **xp, **x1, **loc_X1, **loc_tao;
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
  DMDAVecGetArray(info->da, Y, &loc_tao);

  PetscTime(&v1);
  ComputeLfunction(info, xp, x1, loc_tao, user);                /* compute inner point for x1 */ 
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
  ComputeLfunction(info, loc_X1, x, loc_tao, user); 
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
  DMDAVecRestoreArray(info->da, loc_xp, &xp); 
  DMDAVecRestoreArray(info->da, loc_x1, &loc_X1); 
  DMDAVecRestoreArray(info->da, Y, &loc_tao);
  VecDestroy(&X1); 
  VecDestroy(&loc_xp); 
  VecDestroy(&loc_x1); 
  
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
  user->H = 1;
  user->stime = 0.1 / (user->Nx - 1);                /* for CFL number be small */
  user->rho0 = 1.0; 
  user->maxtime = 200;
  user->W = 0.125;
  user->Re = 500.0;
  user->rate = 99.0;
  user->lidvelocity = 0.1; 
  user->F = 8 * user->lidvelocity * user->rho0 * (100.0/user->Re) / pow(user->H,2);
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
