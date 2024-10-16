#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"

#define LSB 0
#define MSB 1

// From writeVrt.c in strack
//void writeSingleVRT(int32_t nR, int32_t nA, dictNode *metaData, char *vrtFile, char *bandFiles[], char *bandNames[],
//                    GDALDataType dataTypes[], char *byteSwapOption, int32_t nBands);

void writeCullVrt(CullParams *cullPar, int32_t *byteOrder)
{
    char MTVrt[2048], *tmp, *byteSwapOption;
    // Band stuff for MT files
    GDALDataType dataTypesM[] = {GDT_Byte};
    fprintf(stderr, "writing cull vrt...\n");
    sprintf(MTVrt, "%s.vrt", cullPar->outFileT);
    char *bandFilesM[] = {cullPar->outFileT};
    char *bandNamesM[] = {"MatchType"};
    
    fprintf(stderr, "%s\n", get_value(cullPar->metaData, "ByteOrder"));
    // Set byte order option
    if (strstr(get_value(cullPar->metaData, "ByteOrder"), "MSB"))
    {
        byteSwapOption = "BYTEORDER=MSB";
        *byteOrder = MSB;
    }
    else
    {
        byteSwapOption = "BYTEORDER=LSB";
        *byteOrder = LSB;
    }
    // Band info
    GDALDataType dataTypes[] = {GDT_Float32, GDT_Float32, GDT_Float32, GDT_Float32, GDT_Float32};
    char *bandFiles[] = {cullPar->outFileR, cullPar->outFileA, cullPar->outFileSR, cullPar->outFileSA, cullPar->outFileC};
    char *bandNames[] = {"RangeOffsets", "AzimuthOffsets", "RangeSigma", "AzimuthSigma", "Correlation"};
    // Write VRT
    writeSingleVRT(cullPar->nR, cullPar->nA, cullPar->metaData, cullPar->outFileVRT, bandFiles, bandNames, dataTypes, byteSwapOption, -2.0e9, 5);
    // Write Match VRT
    writeSingleVRT(cullPar->nR, cullPar->nA, cullPar->metaDataMT, MTVrt, bandFilesM, bandNamesM, dataTypesM, NULL, DONOTINCLUDENODATA, 1);
}

static size_t fwriteOptionalBS(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags, int32_t byteOrder)
{
    if (byteOrder == LSB)
        return fwrite(ptr, nitems, size, fp);
    else
        return fwriteBS(ptr, nitems, size, fp, flags);
}

void writeCullData(CullParams *cullPar)
{
    FILE *fp;
    int32_t nr, na, byteOrder, i;
    size_t nSamples;
    // Write the vrt
    fprintf(stderr, "writing cull data...\n");
    writeCullVrt(cullPar, &byteOrder);
    nr = cullPar->nR;
    na = cullPar->nA;
    char *files[] = {cullPar->outFileR, cullPar->outFileA, cullPar->outFileSR,
                     cullPar->outFileSA, cullPar->outFileC, cullPar->outFileT};
    void *buffers[] = {cullPar->offRS[0], cullPar->offAS[0], cullPar->sigmaR[0],
                       cullPar->sigmaA[0], cullPar->corr[0], cullPar->type[0]};
    size_t sizes[] = {sizeof(float), sizeof(float), sizeof(float),
                      sizeof(float), sizeof(float), sizeof(char)};
    int32_t flags[] = {FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, BYTEFLAG};
    nSamples = (cullPar->nA) * cullPar->nR;
    for (i = 0; i < 6; i++)
    {
        fprintf(stderr, "Writing %s\n", files[i]);
        fp = fopen(files[i], "w");
        fwriteOptionalBS(buffers[i], sizes[i], nSamples, fp, flags[i], byteOrder);
        fclose(fp);
    }
 
}
