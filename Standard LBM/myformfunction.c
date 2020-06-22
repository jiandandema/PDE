#include <petscdmda.h>
#include "def.h"




PetscErrorCode Collision_Step(DMDALocalInfo *info, Field **x, Field **f, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  PetscScalar    sumrho, sumu, sumv; 


  PetscFunctionBegin;


  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym; 


  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {
        sumrho = 0.0; sumu = 0.0; sumv = 0.0;
        for (k = 0; k < user->Q; k++) {
          sumrho += x[j][i].f[k];
          sumu   += user->eu[k] * x[j][i].f[k];
          sumv   += user->ev[k] * x[j][i].f[k];
        }
        sumu = sumu / sumrho;
        sumv = sumv / sumrho;
        for (k = 0; k < user->Q; k++) {
          f[j][i].f[k] = x[j][i].f[k] - (CalFeq(k, sumrho ,sumu, sumv, user)- x[j][i].f[k]) * / user->tau_star;
        }
      }
    }
  }

  PetscFunctionReturn(0); 
}


PetscErrorCode Propagration_Step(DMDALocalInfo *info, Field **f, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 


  PetscFunctionBegin;


  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym; 


  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {
        for (k = 0; k < user->Q; k++) {
          f[j + user->ev[k]][i + user->eu[k]].f[k] = f[j][i].f[k];
        }
      }
    }
  }

  PetscFunctionReturn(0); 
}


PetscErrorCode ComputeLBoundary(DMDALocalInfo *info, Field **x1, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  PetscScalar    sumrho, sumu, sumv, rho, u, v, pu, pv, prho; 

  PetscFunctionBegin; 

  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym; 

  for (j = yints; j < yinte; j++) {

    if (j == 0)
      for (i = xints; i < xinte; i++) {
        if (i >= 0 && i <= info->mx - 1) {                /* boundary points */
          rho = 0.0;   u = 0.0;   v = 0.0; 
          prho = 0.0;  pu = 0.0;  pv = 0.0;                /* rho, U, V */ 
          for (k = 0; k < user->Q; k++) {
            rho  += x1[j+1][i].f[k]; 
            u    += user->eu[k] * x1[j+1][i].f[k]; 
            v    += user->ev[k] * x1[j+1][i].f[k]; 
            prho += x1[j+2][i].f[k]; 
            pu   += user->eu[k] * x1[j+2][i].f[k]; 
            pv   += user->ev[k] * x1[j+2][i].f[k]; 
          }
          u = u / rho; 
          v = v / rho; 
          pu = pu / prho; 
          pv = pv / prho; 
          sumu = 0.0; 
          sumv = 0.0;                /* the initialguess of u / v in boundary point with boundary condition */
          sumrho = 2 * rho - prho;                /* the initialguess of rho in boundary point */
          for (k = 0; k < user->Q; k++) {
#if highorderext == 1
            x1[j][i].f[k] =  CalFeq(k, sumrho, sumu, sumv, user) + 2 * x1[j+1][i].f[k] - 2 * CalFeq(k, rho, u, v, user)
                             - x1[j+2][i].f[k] + CalFeq(k, prho, pu, pv, user);                /* second order accuracy NEM */
#else          
            x1[j][i].f[k] =  CalFeq(k, sumrho, sumu, sumv, user) + x1[j+1][i].f[k] -  CalFeq(k, rho, u, v, user);
#endif
          }
        }
      }

    if (j == info->my - 1)
      for (i = xints; i < xinte; i++) {
        if (i >= 0 && i <= info->mx - 1) {
          rho = 0.0;   u = 0.0;   v = 0.0; 
          prho = 0.0;  pu = 0.0;  pv = 0.0;                /* rho, U, V */  
          for (k = 0; k < user->Q; k++) {
            rho  += x1[j-1][i].f[k]; 
            u    += user->eu[k] * x1[j - 1][i].f[k]; 
            v    += user->ev[k] * x1[j - 1][i].f[k]; 
            prho += x1[j-2][i].f[k]; 
            pu   += user->eu[k] * x1[j - 2][i].f[k]; 
            pv   += user->ev[k] * x1[j - 2][i].f[k]; 
          }
          u = u / rho; 
          v = v / rho; 
          pu = pu / prho; 
          pv = pv / prho; 
          sumu = user->lidvelocity; 
          sumv = 0.0; 
          sumrho = 2 * rho - prho; 
          for(k = 0; k<user->Q; k++){
#if highorderext == 1
            x1[j][i].f[k] = CalFeq(k, sumrho, sumu, sumv, user) + 2.0 * (x1[j-1][i].f[k] - CalFeq(k, rho, u, v, user))
                            - x1[j-2][i].f[k] + CalFeq(k, prho, pu, pv, user); 
#else          
            x1[j][i].f[k] =  CalFeq(k, sumrho, sumu, sumv, user) + x1[j+1][i].f[k] -  CalFeq(k, rho, u, v, user);
#endif
          }
        }
      }
      
    }

  PetscFunctionReturn(0); 
}
