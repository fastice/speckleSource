#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"
#include "mosaicSource/common/common.h"
#include "stdlib.h"

static float *sWeights(int32_t sA)
{
	int32_t m;
	float w2;  /* half width as float */
	float *wA; /* Weights */
	w2 = (float)sA / 2.0 - 0.01;
	wA = (float *)malloc(sizeof(float) * sA + 1);
	wA += sA/2;
	for (m = -sA / 2; m <= sA / 2; m++)
	{
		if (m < -w2 || m > w2)
			wA[m] = 0.5;
		else
			wA[m] = 1.0;
		// fprintf(stderr, " weights %d %f sA %i\n", m, wA[m], sA);
	}
	return wA;
}

void cullSmooth(CullParams *cullPar)
{
	int32_t i, j;
	int32_t i1, j1, i2, j2;
	int32_t ii, jj;
	int32_t iw, jw, jwStart, iwStart;
	float nlooks[4];
	float *wA, *wR;
	double w;
	double ngood;
	double meanA, meanR;
	float nLooksEff;
	nlooks[0] = 0.0;
	nlooks[1] = 1.0;
	nlooks[2] = 0.25;
	nlooks[3] = 0.015;

	wA = sWeights(cullPar->sA);
	wR = sWeights(cullPar->sR);
	fprintf(stderr, "\n**** smoothing data sr=%i  sa=%i ****\n\n", cullPar->sR, cullPar->sA);
	/* Azimuth loop */
	for (i = 0; i < cullPar->nA; i++)
	{
		/*
		   Compute averaging window
		*/
		i1 = max(0, i - cullPar->sA / 2);
		iwStart = i1 - i;
		i2 = min(cullPar->nA - 1, i + cullPar->sA / 2);
		/* Range loop */
		for (j = 0; j < cullPar->nR; j++)
		{
			j1 = max(0, j - cullPar->sR / 2);
			j2 = min(cullPar->nR - 1, j + cullPar->sR / 2);
			jwStart = j1 - j;
			ngood = 0;
			nLooksEff = 0.0;
			meanA = 0.0;
			meanR = 0.0;
			/* Smoothing window */
			for (ii = i1, iw = iwStart; ii <= i2; ii++, iw++)
			{
				for (jj = j1, jw = jwStart; jj <= j2; jj++, jw++)
				{
					if(jj >= cullPar->nR || ii >= cullPar->nA || jj < 0 || ii <0) fprintf(stderr, "%i %i\n",ii,jj);
					/* Valid data, so add */
					if (cullPar->offR[ii][jj] > (1 - LARGEINT) && cullPar->offA[ii][jj] > (1 - LARGEINT))
					{
						w = wA[iw] * wR[jw];								   /* Determine the combined weight, which is applied to data,ngood, nlooksEff */
						ngood += w;											   /* Update count */
						nLooksEff += w * nlooks[(int)(cullPar->type[ii][jj])]; /* update effect looks based on match type */
						meanA += w * cullPar->offA[ii][jj];					   /* local mean azimuth */
						meanR += w * cullPar->offR[ii][jj];					   /* local mean azimuth */
					}														   /* Endif cullPar.. */
				}															   /* End jj */
			}																   /* End ii */
			/* Added exception for sR=1 or sA=1
			   To save must have, a) ngood > 1 or special case with sr or sa=1, b) valid sigmaA c) valid data for that point
			 */
			if (((ngood > 1) || (cullPar->sA == 1) || (cullPar->sR == 1)) && cullPar->offR[i][j] > (1 - LARGEINT) && cullPar->sigmaA[i][j] > (1 - LARGEINT))
			{
				nLooksEff = max(nLooksEff, 1.);
				cullPar->sigmaA[i][j] /= sqrt((double)nLooksEff); /* reduce local sigma by sqrt(effective looks) */
				cullPar->sigmaR[i][j] /= sqrt((double)nLooksEff);
				cullPar->offRS[i][j] = meanR / (float)ngood; /* Complete mean */
				cullPar->offAS[i][j] = meanA / (float)ngood;
			}
			else
			{
				cullPar->offRS[i][j] = (float)-LARGEINT;
				cullPar->offAS[i][j] = (float)-LARGEINT;
			}
		} /* End for j */
	}	  /* End for i */
}
