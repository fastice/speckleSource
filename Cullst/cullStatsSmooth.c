#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "cRecipes/nrutil.h"
#include <omp.h>

typedef struct {
    float x;
    float y;
} intDataSS;

static void planeCoeffsSS(void *x, int32_t i, double *afunc, int32_t ma)
{
    intDataSS *xy = (intDataSS *)x;
    afunc[1] = xy[0].x;
    afunc[2] = (double)(xy[i].x);
    afunc[3] = (double)(xy[i].y);
}

/* Returns pointer shifted by s/2 for centered indexing; *base_out receives
   the original malloc'd pointer for freeing. */
static float *sWeightsSS(int32_t s, float **base_out)
{
    int32_t m;
    float w2 = (float)s / 2.0f - 0.01f;
    float *base = (float *)malloc(sizeof(float) * (s + 1));
    float *wgt  = base + s / 2;
    for (m = -s / 2; m <= s / 2; m++)
        wgt[m] = (m < -w2 || m > w2) ? 0.5f : 1.0f;
    if (base_out) *base_out = base;
    return wgt;
}

void cullStatsSmooth(CullParams *cullPar)
{
    cullPar->bR = 9;
    cullPar->bA = 9;
    int32_t npts = 1 + (cullPar->bR + 1) * (cullPar->bA + 1);
    int32_t ma = 3;
    float nlooks[4] = {0.0f, 1.0f, 0.25f, 0.015f};

    float *wAbase, *wRbase;
    float *wA = sWeightsSS(cullPar->sA, &wAbase);
    float *wR = sWeightsSS(cullPar->sR, &wRbase);

    fprintf(stderr, "Cull Stats+Smooth\n");
    fprintf(stderr, "\n**** smoothing data sr=%i  sa=%i ****\n\n", cullPar->sR, cullPar->sA);

#pragma omp parallel
    {
        int32_t i, j, i1, i2, j1, j2, ii, jj, k;
        uint32_t ngood;
        double varA, varR, meanA, meanR;
        double *listR  = (double *)malloc(sizeof(double) * npts);
        double *listA  = (double *)malloc(sizeof(double) * npts);
        double *dTmp   = (double *)malloc(sizeof(double) * npts);
        intDataSS *data = (intDataSS *)malloc(sizeof(intDataSS) * npts);
        float oR, oA;
        double aR[5], aA[5];
        double chisq;
        double **u = dmatrix(1, npts + 1, 1, ma);
        double **v = dmatrix(1, ma, 1, ma);
        double  *w = dvector(1, ma);
        int32_t wi1, wi2, wj1, wj2, iw, jw, iwStart, jwStart;
        double ngoodS, wt, meanAS, meanRS;
        float nLooksEff;

        for (k = 0; k < npts; k++) dTmp[k] = 1.0;

#pragma omp for schedule(dynamic)
        for (i = 0; i < cullPar->nA; i++)
        {
            /* stats box bounds for this row */
            i1 = max(0, i - cullPar->bA / 2);
            i2 = min(cullPar->nA - 1, i + cullPar->bA / 2);
            /* smooth box bounds for this row */
            wi1 = max(0, i - cullPar->sA / 2);
            iwStart = wi1 - i;
            wi2 = min(cullPar->nA - 1, i + cullPar->sA / 2);

            for (j = 0; j < cullPar->nR; j++)
            {
                /* stats box bounds for this column */
                j1 = max(0, j - cullPar->bR / 2);
                j2 = min(cullPar->nR - 1, j + cullPar->bR / 2);

                /* --- STATS: fit plane, compute local variance --- */
                ngood = 0;
                meanA = 0.0;
                meanR = 0.0;
                aR[1] = aR[2] = aR[3] = 0.0;
                aA[1] = aA[2] = aA[3] = 0.0;
                for (ii = i1; ii <= i2; ii++)
                    for (jj = j1; jj <= j2; jj++)
                        if (cullPar->offR[ii][jj] > (1 - LARGEINT) &&
                            cullPar->offA[ii][jj] > (1 - LARGEINT))
                        {
                            ngood++;
                            listR[ngood] = cullPar->offR[ii][jj];
                            listA[ngood] = cullPar->offA[ii][jj];
                            meanA += cullPar->offA[ii][jj];
                            meanR += cullPar->offR[ii][jj];
                            data[ngood].x = (float)ii;
                            data[ngood].y = (float)jj;
                        }
                if (ngood > 8)
                {
                    meanA /= (double)ngood;
                    meanR /= (double)ngood;
                    data[0].x = meanR;
                    svdfit((void *)data, listR, dTmp, ngood, aR, ma, u, v, w,
                           &chisq, &planeCoeffsSS);
                    data[0].x = meanA;
                    svdfit((void *)data, listA, dTmp, ngood, aA, ma, u, v, w,
                           &chisq, &planeCoeffsSS);
                    aR[1] *= meanR;
                    aA[1] *= meanA;
                }
                ngood = 0;
                meanA = 0.0;
                meanR = 0.0;
                varR  = 0.0;
                varA  = 0.0;
                for (ii = i1; ii <= i2; ii++)
                    for (jj = j1; jj <= j2; jj++)
                        if (cullPar->offR[ii][jj] > (1 - LARGEINT) &&
                            cullPar->offA[ii][jj] > (1 - LARGEINT))
                        {
                            ngood++;
                            oR = cullPar->offR[ii][jj];
                            oA = cullPar->offA[ii][jj];
                            oR -= (aR[1] + aR[2] * (float)ii + aR[3] * (float)jj);
                            oA -= (aA[1] + aA[2] * (float)ii + aA[3] * (float)jj);
                            meanA += oA;
                            meanR += oR;
                            varA  += oA * oA;
                            varR  += oR * oR;
                        }
                if (ngood > 8 && cullPar->offR[i][j] > (1 - LARGEINT))
                {
                    cullPar->sigmaA[i][j] = varA / (float)(ngood - 3);
                    cullPar->sigmaR[i][j] = varR / (float)(ngood - 3);
                    cullPar->sigmaA[i][j] -= meanA * meanA /
                                             (float)((ngood - 3) * ngood);
                    cullPar->sigmaR[i][j] -= meanR * meanR /
                                             (float)((ngood - 3) * ngood);
                    cullPar->sigmaA[i][j] = max(0.00001, cullPar->sigmaA[i][j]);
                    cullPar->sigmaR[i][j] = max(0.00001, cullPar->sigmaR[i][j]);
                    cullPar->sigmaA[i][j] = sqrt((double)cullPar->sigmaA[i][j]);
                    cullPar->sigmaR[i][j] = sqrt((double)cullPar->sigmaR[i][j]);
                }
                else
                {
                    cullPar->sigmaA[i][j] = (float)-LARGEINT;
                    cullPar->sigmaR[i][j] = (float)-LARGEINT;
                }

                /* --- SMOOTH: weighted window average --- */
                wj1 = max(0, j - cullPar->sR / 2);
                wj2 = min(cullPar->nR - 1, j + cullPar->sR / 2);
                jwStart = wj1 - j;
                ngoodS    = 0.0;
                nLooksEff = 0.0f;
                meanAS    = 0.0;
                meanRS    = 0.0;
                for (ii = wi1, iw = iwStart; ii <= wi2; ii++, iw++)
                    for (jj = wj1, jw = jwStart; jj <= wj2; jj++, jw++)
                        if (cullPar->offR[ii][jj] > (1 - LARGEINT) &&
                            cullPar->offA[ii][jj] > (1 - LARGEINT))
                        {
                            wt = wA[iw] * wR[jw];
                            ngoodS    += wt;
                            nLooksEff += wt * nlooks[(int)(cullPar->type[ii][jj])];
                            meanAS    += wt * cullPar->offA[ii][jj];
                            meanRS    += wt * cullPar->offR[ii][jj];
                        }
                if (((ngoodS > 1.0) || (cullPar->sA == 1) || (cullPar->sR == 1)) &&
                    cullPar->offR[i][j] > (1 - LARGEINT) &&
                    cullPar->sigmaA[i][j] > (1 - LARGEINT))
                {
                    nLooksEff = max(nLooksEff, 1.0f);
                    cullPar->sigmaA[i][j] /= sqrt((double)nLooksEff);
                    cullPar->sigmaR[i][j] /= sqrt((double)nLooksEff);
                    cullPar->offRS[i][j] = (float)(meanRS / ngoodS);
                    cullPar->offAS[i][j] = (float)(meanAS / ngoodS);
                }
                else
                {
                    cullPar->offRS[i][j] = (float)-LARGEINT;
                    cullPar->offAS[i][j] = (float)-LARGEINT;
                }
            } /* end for j */
        }     /* end for i */

        free(listR); free(listA); free(dTmp); free(data);
        free_dmatrix(u, 1, npts + 1, 1, ma);
        free_dmatrix(v, 1, ma, 1, ma);
        free_dvector(w, 1, ma);
    } /* end parallel */

    free(wAbase);
    free(wRbase);
}
