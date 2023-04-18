#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"
#include "stdlib.h"
#include "ers1/ers1Code/ers1.h"

static float selectForCull(uint32_t k, uint32_t n, float arr[]);

void cullSTData(CullParams *cullPar)
{
    int32_t i, j;
    int32_t i1, j1, i2, j2;
    int32_t ii, jj;
    uint32_t ngood;
    int32_t nCulled;
    uint32_t midIndex;
    double diffR, diffA;
    float medianR, medianA;
    float *listR, *listA;
    /*    cullPar->bR=9;
      cullPar->bA=9; moved to readargs */
    listR = (float *)
        malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));
    listA = (float *)
        malloc(sizeof(float) * (100 + (cullPar->bR + 1) * (cullPar->bA + 1)));

    nCulled = 0;
    ngood = 0;
    for (i = 0; i < cullPar->nA; i++)
    {
        for (j = 0; j < cullPar->nR; j++)
        {
            if (cullPar->offR[i][j] > (1 - LARGEINT) &&
                cullPar->offA[i][j] > (1 - LARGEINT))
                ngood++;
            else
                nCulled++;
        }
    }

    fprintf(stderr, "n initial %f\n", (double)nCulled / (double)(nCulled + ngood));
    ngood = 0;
    nCulled = 0;
    for (i = 0; i < cullPar->nA; i++)
    {
        i1 = max(0, i - cullPar->bA / 2);
        i2 = min(cullPar->nA - 1, i + cullPar->bA / 2);
        for (j = 0; j < cullPar->nR; j++)
        {
            j1 = max(0, j - cullPar->bR / 2);
            j2 = min(cullPar->nR - 1, j + cullPar->bR / 2);

            ngood = 0;
            for (ii = i1; ii <= i2; ii++)
            {
                for (jj = j1; jj <= j2; jj++)
                {
                    if (cullPar->offR[ii][jj] > (1 - LARGEINT) &&
                        cullPar->offA[ii][jj] > (1 - LARGEINT))
                    {
                        ngood++;
                        listR[ngood] = cullPar->offR[ii][jj];
                        listA[ngood] = cullPar->offA[ii][jj];
                    } /* Endif cullPar.. */
                }     /* End jj */
            }         /* End ii */
            if (ngood > cullPar->nGood)
            {
                midIndex = ngood / 2;
                medianR = selectForCull(midIndex, ngood, listR);
                medianA = selectForCull(midIndex, ngood, listA);
                diffA = (double)(cullPar->offA[i][j] - medianA);
                diffR = (double)(cullPar->offR[i][j] - medianR);
                if (fabs(diffA) > cullPar->maxA || fabs(diffR) > cullPar->maxR)
                {
                    cullPar->offR[i][j] = (float)-LARGEINT;
                    cullPar->offA[i][j] = (float)-LARGEINT;
                    nCulled++;
                }
            }
            else
            {
                cullPar->offR[i][j] = (float)-LARGEINT;
                cullPar->offA[i][j] = (float)-LARGEINT;
                nCulled++;
            }
        } /* End for j */
    }     /* End for i */
    fprintf(stderr, "nCulled %i %f\n", nCulled,
            (double)nCulled / (double)(cullPar->nR * cullPar->nA));
}

#define SWAP(a, b) \
    temp = (a);    \
    (a) = (b);     \
    (b) = temp;

float selectForCull(uint32_t k, uint32_t n, float arr[])
{
    uint32_t i, ir, j, l, mid;
    float a, temp;

    l = 1;
    ir = n;
    for (;;)
    {
        if (ir <= l + 1)
        {
            if (ir == l + 1 && arr[ir] < arr[l])
            {
                SWAP(arr[l], arr[ir])
            }
            return arr[k];
        }
        else
        {
            mid = (l + ir) >> 1;
            SWAP(arr[mid], arr[l + 1])
            if (arr[l] > arr[ir])
            {
                SWAP(arr[l], arr[ir])
            }
            if (arr[l + 1] > arr[ir])
            {
                SWAP(arr[l + 1], arr[ir])
            }
            if (arr[l] > arr[l + 1])
            {
                SWAP(arr[l], arr[l + 1])
            }
            i = l + 1;
            j = ir;
            a = arr[l + 1];
            for (;;)
            {
                do
                    i++;
                while (arr[i] < a);
                do
                    j--;
                while (arr[j] > a);
                if (j < i)
                    break;
                SWAP(arr[i], arr[j])
            }
            arr[l + 1] = arr[j];
            arr[j] = a;
            if (j >= k)
                ir = j - 1;
            if (j <= k)
                l = i;
        }
    }
}
#undef SWAP
