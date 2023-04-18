#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"
#include <stdlib.h>

static char **mallocByteMat(int32_t nA, int32_t nR);
static float **mallocFloatMat(int32_t nA, int32_t nR);

void loadCullData(CullLSParams *cullPar)
{
    FILE *fp;
    char line[1024];
    int32_t lineCount, eod;
    int32_t r0, a0, nr, na;
    int32_t deltaA, deltaR;
    /*
         Read .dat file
    */
    fp = openInputFile(cullPar->inFileD);
    lineCount = getDataString(fp, lineCount, line, &eod);
    if (sscanf(line, "%i%i%i%i%i%i", &r0, &a0, &nr, &na, &deltaR, &deltaA) != 6)
        error("%s  %i of %s\n%s",
              "readOffsets -- Missing image parameters at line:",
              lineCount, cullPar->inFileD, line);
    close(fp);
    cullPar->nR = nr;
    cullPar->nA = na;
    cullPar->rStart = r0;
    cullPar->aStart = a0;
    cullPar->deltaA = deltaA;
    cullPar->deltaR = deltaR;

    /*
       Mallocldatam
     */
    cullPar->offR = mallocFloatMat(na, nr);
    cullPar->offA = mallocFloatMat(na, nr);
    cullPar->offRS = mallocFloatMat(na, nr);
    cullPar->offAS = mallocFloatMat(na, nr);
    cullPar->corr = mallocFloatMat(na, nr);
    cullPar->sigmaR = mallocFloatMat(na, nr);
    cullPar->sigmaA = mallocFloatMat(na, nr);
    cullPar->type = mallocByteMat(na, nr);
    /*
       Load data
    */
    fp = openInputFile(cullPar->inFileR);
    freadBS(cullPar->offR[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fp = openInputFile(cullPar->inFileA);
    freadBS(cullPar->offA[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fp = openInputFile(cullPar->inFileC);
    freadBS(cullPar->corr[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fp = openInputFile(cullPar->inFileT);
    freadBS(cullPar->type[0], sizeof(char), nr * na, fp, BYTEFLAG);
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
