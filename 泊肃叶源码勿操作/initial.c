#include <petscdmda.h>
#include "def.h"


PetscErrorCode Getanalysis(DMDALocalInfo *info, AppCtx *user)
{
  PetscInt       i, j, k;
  Field          **x;
  PetscScalar    xj;

  DMDAVecGetArray(info->da, user->func, &x);

  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      for (k = 1; k < user->Q; k++) {
        x[j][i].f[k] = 0.0;
      }
      xj = 1.0 * ((PetscScalar)(j)) / (info->my - 1);
      x[j][i].f[0] = user->lidvelocity * 4 * xj * (1 - xj);                /* real solution:rho = user->rho0; u = user->lidvelocity * 4 * xj * (1-xj); v = 0.0; */
    }
  }

  DMDAVecRestoreArray(info->da, user->func, &x);

  return 0;
}


PetscScalar Error(DMDALocalInfo *info, Vec X, AppCtx *user)
{
  PetscInt       i, j, k;
  Field          **px, **x, **error, **Uv;
  Vec            ex, uv;
  PetscScalar    sumu[2], err, sumrho, u, v, pu, pv;

  VecDuplicate(X, &ex);
  VecDuplicate(X, &uv);

  DMDAVecGetArray(info->da, X, &x);
  DMDAVecGetArray(info->da, user->func, &px);
  DMDAVecGetArray(info->da, ex, &error);
  DMDAVecGetArray(info->da, uv, &Uv);

  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      for (k = 2; k < user->Q; k++) {
        error[j][i].f[k] = 0.0;
        Uv[j][i].f[k] = 0.0;
      }
      
      sumrho = 0; 
      u = 0; 
      v = 0;
      for(k = 0; k < user->Q; k++){
        sumrho += x[j][i].f[k];
        u += user->eu[k] * x[j][i].f[k];
      }
      u  = u / sumrho; 
      pu = px[j][i].f[0];                /* px[j][i].f[0] refer to Getanalysis */
      v  = 0; 
      pv = 0;

      error[j][i].f[0] = u - pu;
      Uv[j][i].f[0] = pu;
      error[j][i].f[1] = v - pv;
      Uv[j][i].f[1] = pv;
    }
  }

  DMDAVecRestoreArray(info->da, X, &x);
  DMDAVecRestoreArray(info->da, user->func, &px);
  DMDAVecRestoreArray(info->da, ex, &error);
  DMDAVecRestoreArray(info->da, uv, &Uv);

  sumu[0] = 0.0; 
  sumu[1] = 0.0;
  VecNorm(uv, NORM_2, sumu);
  VecNorm(ex, NORM_2, sumu+1);
  err = (sumu[1] / (sumu[0] + 1e-10));
  user->ee[i-1][0] = err;
  PetscPrintf(PETSC_COMM_WORLD, "current L^2 error = %g  %g  %g\n", err, sumu[0], sumu[1]);

  VecDestroy(&ex);
  VecDestroy(&uv);
  return err;
}


PetscErrorCode DataSaveASCII(DMDALocalInfo *info, Vec X, char *filename, char *remakefilename)
{
  PetscViewer    dataviewer;                /* http://www.mcs.anl.gov/petsc/petsc-3.2-p7/docs/manualpages/Viewer/PetscViewer.html */
  double         line_solution, mid_solution[360000];
  int            m, n = 360000;
  FILE           *fpf, *fpt;

  PetscFunctionBegin;                /* First executable line of each PETSc function used for error handling */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &dataviewer);                 /* Opens an ASCII file as a PetscViewer */
  PetscViewerPushFormat(dataviewer, PETSC_VIEWER_ASCII_SYMMODU);                 /* http://www.mcs.anl.gov/petsc/petsc-3.2-p7/docs/manualpages/Viewer/PetscViewerSetFormat.html */
  VecView(X, dataviewer);                 /* Views a vector object */
  PetscViewerDestroy(&dataviewer);                 /* Destroys a PetscViewer */

  fpf = fopen(filename, "r");
  for (m = 0; m < 2; m++) 
  	fscanf(fpf, "%*[^\n]%*c"); 
  for (m = 0; m < n; m++) {
    fscanf(fpf, "%lE", &line_solution);
    mid_solution[m] = line_solution;
  }
  fclose(fpf);

  fpt = fopen(remakefilename, "w");
  for (m = 0; m < n; m++)
    fprintf(fpt, "%.16e\n", mid_solution[m]);
  fclose(fpt);

  PetscFunctionReturn(0);                /* Last executable line of each PETSc function used for error handling. Replaces return() */
}
/*******************************************************************************************
* YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
********************************************************************************************/
PetscErrorCode DataSaveBin(Vec x, char *filename)
{
  PetscViewer    dataviewer;

  PetscFunctionBegin;

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &dataviewer);                 /* Opens a file for binary input/output */
  VecView(x, dataviewer); 
  PetscViewerDestroy(&dataviewer); 

  PetscFunctionReturn(0);
}
