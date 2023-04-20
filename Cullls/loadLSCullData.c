#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include <stdlib.h>
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "cullls.h"
static char **mallocByteMat(int32_t nA, int32_t nR);
static float **mallocFloatMat(int32_t nA, int32_t nR);

void loadLSCullData(CullLSParams *cullPar)
{
	int32_t nx, ny;
	/*
	  Read offsets
	 */
	readLSOffsets(&(cullPar->fitDat), &(cullPar->matches), TRUE, NULL);
	/*
malloc space for sigma buffers
	 */
	nx = cullPar->matches.nx;
	ny = cullPar->matches.ny;
	cullPar->matches.sigmaX = mallocFloatMat(ny, nx);
	cullPar->matches.sigmaY = mallocFloatMat(ny, nx);
	/*
	  Malloc other data ??? verify we need these
	*/
	cullPar->offXC = mallocFloatMat(ny, nx);
	cullPar->offYC = mallocFloatMat(ny, nx);
	cullPar->XS = mallocFloatMat(ny, nx);
	cullPar->YS = mallocFloatMat(ny, nx);
	/*  cullPar->corr = mallocFloatMat(ny,nx);
		cullPar->type = mallocByteMat(ny,nx);*/
	/*
	  Load data
	*/
}

static char **mallocByteMat(int32_t nA, int32_t nR)
{
	char *tmp, **tmp1;

	int32_t i;
	tmp = malloc(nR * nA * sizeof(char));
	tmp1 = (char **)malloc(nA * sizeof(char *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}

static float **mallocFloatMat(int32_t nA, int32_t nR)
{
	float *tmp, **tmp1;
	int32_t i;
	size_t size;
	size = (size_t)nR * (size_t)nA * sizeof(float);
	tmp = malloc(size);
	tmp1 = (float **)malloc(nA * sizeof(float *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}
