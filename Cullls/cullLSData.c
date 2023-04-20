#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "cullls.h"
/*#include "ers1/ers1Code/ers1.h"*/
static float selectForCull(uint32_t k, uint32_t n, float arr[]);

void cullLSData(CullLSParams *cullPar)
{
	int32_t i, j;
	int32_t i1, j1, i2, j2;
	int32_t ii, jj;
	uint32_t nGoodPts, ngood;
	int32_t nCulled, nWithData, nTotal;
	uint32_t midIndex;
	double diffX, diffY;
	float medianX, medianY;
	float *listX, *listY;
	double maxX, maxY;
	int32_t slowFlag;
	int32_t nx, ny;
	fprintf(stderr, "Entering cullLSData\n");
	listX = (float *)
		malloc(sizeof(float) * (1000 + (cullPar->bR + 1) * (cullPar->bA + 1)));
	listY = (float *)
		malloc(sizeof(float) * (100 + (cullPar->bR + 1) * (cullPar->bA + 1)));
	nCulled = 0;
	nGoodPts = 0;
	nx = cullPar->matches.nx;
	ny = cullPar->matches.ny;
	nWithData = 0;
	nTotal = 0;
	for (i = 0; i < ny; i++)
	{
		for (j = 0; j < nx; j++)
		{
			if (cullPar->matches.X[i][j] > (1 - LARGEINT) &&
				cullPar->matches.Y[i][j] > (1 - LARGEINT))
				nGoodPts++;
			else
				nCulled++;
			nTotal++;
			if ((cullPar->matches.type[i][j] & NODATATOUSE) == 0)
				nWithData++;
		}
	}
	/*
	   nMatch set in main to -1, after that don't save (i.e, save original stat )
	 */
	if (cullPar->nMatch < 0)
	{
		cullPar->nAttempt = nWithData;
		cullPar->nMatch = nGoodPts;
	}
	fprintf(stderr, "n initial %f %f %i %i %i %i\n", (double)nGoodPts / (double)(nCulled + nGoodPts),
			(float)nGoodPts / (float)nWithData, nTotal, nWithData, ngood, nCulled);

	ngood = 0;
	nCulled = 0;
	nGoodPts = 0;
	for (i = 0; i < ny; i++)
	{
		i1 = max(0, i - cullPar->bA / 2);
		i2 = min(ny - 1, i + cullPar->bA / 2);
		for (j = 0; j < nx; j++)
		{
			j1 = max(0, j - cullPar->bR / 2);
			j2 = min(nx - 1, j + cullPar->bR / 2);
			slowFlag = TRUE;
			ngood = 0;
			for (ii = i1; ii <= i2; ii++)
			{
				for (jj = j1; jj <= j2; jj++)
				{
					if (cullPar->matches.X[ii][jj] > (1 - LARGEINT) &&
						cullPar->matches.Y[ii][jj] > (1 - LARGEINT))
					{
						ngood++;
						listX[ngood] = cullPar->matches.X[ii][jj];
						listY[ngood] = cullPar->matches.Y[ii][jj];
						if (!(cullPar->matches.type[ii][jj] & SLOWVEL))
							slowFlag = FALSE;
					} /* Endif cullPar.. */
				}	  /* End jj */
			}		  /* End ii */
			if (ngood > cullPar->nGood && cullPar->matches.Y[i][j] > (NODATA + 1) && cullPar->matches.Y[i][j] > (NODATA + 1))
			{
				midIndex = ngood / 2;
				medianX = selectForCull(midIndex, ngood, listX);
				medianY = selectForCull(midIndex, ngood, listY);
				diffY = (double)(cullPar->matches.Y[i][j] - medianY);
				diffX = (double)(cullPar->matches.X[i][j] - medianX);

				if (slowFlag == TRUE)
				{
					/* scale for time, based on 40 m/yr max deviation */
					maxX = 40.0 / 365. * cullPar->fitDat.deltaT;
					maxY = maxX;
					/* Scale for pixel size */
					maxX /= cullPar->matches.dx;
					maxY /= cullPar->matches.dy;
					maxX = max(maxX, 0.2); /* Allow at least +/- 0.2 */
					maxY = max(maxY, 0.2);
					/*					fprintf(stderr,"%lf %lf %lf %lf %lf %lf  %lf\n",diffX,diffY,cullPar->fitDat.deltaT,medianX,medianY,maxX,maxY);*/
				}
				else
				{
					maxX = cullPar->maxX;
					maxY = cullPar->maxY;
				}
				if (fabs(diffY) > maxY || fabs(diffX) > maxX)
				{
					cullPar->matches.X[i][j] = (float)-LARGEINT;
					cullPar->matches.Y[i][j] = (float)-LARGEINT;
					nCulled++;
				}
				else
					nGoodPts++;
			}
			else
			{
				cullPar->matches.X[i][j] = (float)-LARGEINT;
				cullPar->matches.Y[i][j] = (float)-LARGEINT;
				nCulled++;
			}
		} /* End for j */
	}	  /* End for i */
	fprintf(stderr, "n initial %f %f %i %i %i %i\n", (double)nGoodPts / (double)(nCulled + nGoodPts),
			(float)nGoodPts / (float)nWithData, nTotal, nWithData, ngood, nCulled);
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
