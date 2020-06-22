#include <petscdmda.h>
#include "def.h"


PetscErrorCode ComputeLfunction(DMDALocalInfo *info, Field **x, Field **f, AppCtx *user)
{
  PetscInt       xints, xinte, yints, yinte, i, j, k; 
  PetscScalar    sumrho, sumu, sumv, F_k; 
  PetscScalar    Dfx, Dfy, dx, dy, epsilon; 

  PetscFunctionBegin;                /* First executable line of each PETSc function, used for error handling. Final line of PETSc functions should be PetscFunctionReturn(0) */ 

  xints = info->xs; xinte = info->xs + info->xm; yints = info->ys; yinte = info->ys + info->ym; 
  dx = (BOUNDARYX==1)?(1.0 / info->mx):(1.0 / (info->mx - 1));                /* ?(DMDA_BOUNDARY_PERIODIC):(DMDA_BOUNDARY_NONE) */
  dy = (BOUNDARYY==1)?(1.0 / info->my):(1.0 / (info->my - 1)); 
  epsilon = user->ep; 

  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      if (j > 0 && j < info->my - 1) {
        sumrho = 0.0; sumu = 0.0; sumv = 0.0;
        for (k = 0; k < user->Q; k++) {
          sumrho += x[j][i].f[k];
          sumu   += user->eu[k] * x[j][i].f[k];
          sumv   += user->ev[k] * x[j][i].f[k];
        }
        sumu += 0.5 * user->F_x * user->stime;
        sumv += 0.5 * user->F_y * user->stime;
        sumu = sumu / sumrho;
        sumv = sumv / sumrho;
        for (k = 0; k < user->Q; k++) {
          Dfx = user->eu[k] * (x[j][i+1].f[k] - x[j][i-1].f[k]) / 2.0 / dx;
          Dfy = user->ev[k] * (x[j+1][i].f[k] - x[j-1][i].f[k]) / 2.0 / dy;
          F_k = (1.0 - user->stime / (2.0 * user->tau)) * user->w[k] * (3 * ((user->eu[k] - sumu) * user->F_x + (user->ev[k] - sumv) * user->F_y) + 9 * (user->eu[k] * sumu + user->ev[k] * sumv) * (user->eu[k] * user->F_x + user->ev[k] * user->F_y));
          f[j][i].f[k] = -epsilon * Dfx - epsilon * Dfy + (CalFeq(k, sumrho ,sumu, sumv, user)- x[j][i].f[k]) / user->tau + F_k;
          if (user->eu[k] >= 0) {
            Dfx = user->eu[k] * (3 * x[j][i].f[k] - 4 * x[j][i-1].f[k] + x[j][i-2].f[k]) / 2.0 / dx;                /* upwind scheme */
          }
          else {
            Dfx = -user->eu[k] * (3 * x[j][i].f[k] - 4 * x[j][i+1].f[k] + x[j][i+2].f[k]) / 2.0 / dx;
          }
          if (user->ev[k] >= 0) {
            if (j == 1)
              Dfy = user->ev[k] * (x[j+2][i].f[k] - 2 * x[j+1][i].f[k] 
                  + 3 * x[j][i].f[k] - 2 * x[j-1][i].f[k]) / 2.0 / dy;                /* upwind scheme *//* interpolation for virtual layer */
            else
              Dfy = user->ev[k] * (3 * x[j][i].f[k] - 4 * x[j-1][i].f[k] + x[j-2][i].f[k]) / 2.0 / dy;
          }
          else {
            if (j == info->my - 2)
              Dfy = -user->ev[k] * (x[j-2][i].f[k] - 2 * x[j-1][i].f[k] 
                  + 3 * x[j][i].f[k] - 2 * x[j+1][i].f[k]) / 2.0 / dy;                /* upwind scheme *//* interpolation for virtual layer */
            else
              Dfy = -user->ev[k] * (3 * x[j][i].f[k] - 4 * x[j+1][i].f[k] + x[j+2][i].f[k]) / 2.0 / dy;
          }
          f[j][i].f[k] += -(1- epsilon) * Dfx - (1-epsilon) * Dfy;
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
          sumu = 0.0; 
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
