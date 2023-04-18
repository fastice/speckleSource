#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"

void writeCullData(CullParams *cullPar)
{
    FILE *fp;
    int32_t nr, na;
    nr = cullPar->nR;
    na = cullPar->nA;
    /*
       Load data
    */
    fp = fopen(cullPar->outFileR, "w");
    fwriteBS(cullPar->offRS[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fclose(fp);
    fp = fopen(cullPar->outFileA, "w");
    fwriteBS(cullPar->offAS[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fclose(fp);

    fp = fopen(cullPar->outFileC, "w");
    fwriteBS(cullPar->corr[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fclose(fp);
    fp = fopen(cullPar->outFileT, "w");
    fwriteBS(cullPar->type[0], sizeof(char), nr * na, fp, BYTEFLAG);
    fclose(fp);

    fp = fopen(cullPar->outFileSR, "w");
    fwriteBS(cullPar->sigmaR[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fclose(fp);

    fp = fopen(cullPar->outFileSA, "w");
    fwriteBS(cullPar->sigmaA[0], sizeof(float), nr * na, fp, FLOAT32FLAG);
    fclose(fp);
}
