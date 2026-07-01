#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include <omp.h>
#define NARROWREGION 4

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
	unsigned char maskVal;

	/*    cullPar->bR=9;
		  cullPar->bA=9; moved to readargs */
	listR = (float *)malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));
	listA = (float *)malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));

	nCulled = 0;
	ngood = 0;
	fprintf(stderr, "maskFlag %d\n", cullPar->maskFlag);
#pragma omp parallel for collapse(2) private(maskVal) reduction(+:ngood,nCulled)
	for (i = 0; i < cullPar->nA; i++)
	{
		for (j = 0; j < cullPar->nR; j++)
		{
			if(cullPar->corr[i][j] < cullPar->corrThresh)
			{
				cullPar->offR[i][j] = -LARGEINT;
				cullPar->offA[i][j] = -LARGEINT;
				cullPar->type[i][j] = 0;
			}
			if (cullPar->maskFlag == TRUE)
			{
				/* Added 3/25/2017 to cull mt3 in regions flagged as narrow */
				maskVal = (unsigned char)cullPar->mask[i][j] & NARROWREGION;
				if ( (maskVal > 0 && cullPar->type[i][j] > 2) )
				{ /* This will get rid of wide windows in narrow regions */
				
					cullPar->offR[i][j] = -LARGEINT;
					cullPar->offA[i][j] = -LARGEINT;
					cullPar->type[i][j] = 0;
				}
				/* added 10/12/2018 for evaluating offsets in tests - keeps only a single match type */
				if ((cullPar->singleMT > 0) && (cullPar->type[i][j] != cullPar->singleMT))
				{
					cullPar->offR[i][j] = -LARGEINT;
					cullPar->offA[i][j] = -LARGEINT;
					cullPar->type[i][j] = 0;
				}
			}
			if (cullPar->offR[i][j] > (1 - LARGEINT) && cullPar->offA[i][j] > (1 - LARGEINT))
				ngood++;
			else
				nCulled++;
		}
	}

	fprintf(stderr, "n initial %f\n", (double)nCulled / (double)(nCulled + ngood));

	/* Snapshot offR/offA so all threads read from a consistent pre-pass state */
	float **offRsnap = mallocImage(cullPar->nA, cullPar->nR);
	float **offAsnap = mallocImage(cullPar->nA, cullPar->nR);
	for (i = 0; i < cullPar->nA; i++)
	{
		memcpy(offRsnap[i], cullPar->offR[i], cullPar->nR * sizeof(float));
		memcpy(offAsnap[i], cullPar->offA[i], cullPar->nR * sizeof(float));
	}

	nCulled = 0;
#pragma omp parallel
	{
		int32_t ti, tj, ti1, ti2, tj1, tj2, tii, tjj;
		uint32_t tngood, tmidIndex;
		double tdiffR, tdiffA;
		float tmedianR, tmedianA;
		int32_t tnCulled = 0;
		float *tlistR = (float *)malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));
		float *tlistA = (float *)malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));

#pragma omp for schedule(dynamic)
		for (ti = 0; ti < cullPar->nA; ti++)
		{
			ti1 = max(0, ti - cullPar->bA / 2);
			ti2 = min(cullPar->nA - 1, ti + cullPar->bA / 2);
			for (tj = 0; tj < cullPar->nR; tj++)
			{
				tj1 = max(0, tj - cullPar->bR / 2);
				tj2 = min(cullPar->nR - 1, tj + cullPar->bR / 2);
				tngood = 0;
				for (tii = ti1; tii <= ti2; tii++)
				{
					for (tjj = tj1; tjj <= tj2; tjj++)
					{
						if (offRsnap[tii][tjj] > (1 - LARGEINT) && offAsnap[tii][tjj] > (1 - LARGEINT))
						{
							tngood++;
							tlistR[tngood] = offRsnap[tii][tjj];
							tlistA[tngood] = offAsnap[tii][tjj];
						}
					}
				}
				if (tngood > (uint32_t)cullPar->nGood)
				{
					tmidIndex = tngood / 2;
					tmedianR = selectForCull(tmidIndex, tngood, tlistR);
					tmedianA = selectForCull(tmidIndex, tngood, tlistA);
					tdiffA = (double)(cullPar->offA[ti][tj] - tmedianA);
					tdiffR = (double)(cullPar->offR[ti][tj] - tmedianR);
					if (fabs(tdiffA) > cullPar->maxA || fabs(tdiffR) > cullPar->maxR)
					{
						cullPar->offR[ti][tj] = (float)-LARGEINT;
						cullPar->offA[ti][tj] = (float)-LARGEINT;
						tnCulled++;
					}
				}
				else
				{
					cullPar->offR[ti][tj] = (float)-LARGEINT;
					cullPar->offA[ti][tj] = (float)-LARGEINT;
					tnCulled++;
				}
			}
		}
		free(tlistR);
		free(tlistA);
#pragma omp atomic
		nCulled += tnCulled;
	} /* end parallel */

	for (i = 0; i < cullPar->nA; i++) { free(offRsnap[i]); free(offAsnap[i]); }
	free(offRsnap);
	free(offAsnap);

	fprintf(stderr, "nCulled %i %f\n", nCulled,
			(double)nCulled / (double)(cullPar->nR * cullPar->nA));
}

#define SWAP(a, b) \
	temp = (a);    \
	(a) = (b);     \
	(b) = temp;

float selectForCull(uint32_t k, uint32_t n, float arr[])
// Sort until middle value found
// k = midindex, n = nGood, arr = list
{
	uint32_t i, ir, j, l, mid, m;
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
			// get middle index
			mid = (l + ir) >> 1;
			if(arr[l+1] > arr[mid])
				SWAP(arr[mid], arr[l + 1]);
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
