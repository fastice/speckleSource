#include <stdio.h>


typedef struct
{
    int32_t rStart; /* Starting location for first estimate */
    int32_t aStart;
    int32_t deltaR; /* Jump deltaR/A pixels between estimates */
    int32_t deltaA;
    int32_t nR; /* Number of estimates in range/azimuth direction (nR*nA tot) */
    int32_t nA;
    int32_t maskFlag;
    int32_t useSim;
    int32_t singleMT;
    float **offR;
    float **offA;
    float **offSimA;
    float **offSimR;
    float **offRS; /* SMoothed data */
    float **offAS;
    float **corr;
    float **sigmaA;
    float **sigmaR;
    char **type;
    char **mask;
    char *outFileR;  /* Range ouput file */
    char *outFileA;  /* Azimuth ouput file */
    char *outFileC;  /* Correlation ouput file */
    char *outFileT;  /* Match type ouput file */
    char *outFileD;  /* Header (.dat) file */
    char *outFileVRT;  /* Header (.vrt) file */
     char *outFileVRTMT;  /* Header (.vrt) file */
    char *outFileSR; /* standard deviations */
    char *outFileSA;
    char *inFileR;    /* Range input file */
    char *inFileA;    /* Azimuth input file */
    char *inFileC;    /* Correlation input file */
    char *inFileT;    /* Match type input file */
    char *inFileD;    /* Header (.dat) file */
    char *inFileVRT;    /* Header (.dat) file */
    char *inFileVRTMT;    /* Header (.dat) file */
    char *inFileSimR; /* Range input file */
    char *inFileSimA; /* Azimuth input file */
    char *inFileSimD; /* Header (.dat) file */
    char *inFileSimVRT; /* Header (.dat) file */
    char *geo1;
    char *geo2;
    int32_t bR;       /* Size of window used for checking */
    int32_t bA;
    int32_t nGood;
    float maxA, maxR;
    int32_t sR; /* smoothing window size */
    int32_t sA;
    int32_t islandThresh;
    int32_t ignoreOffsets;
    dictNode *metaData;
    dictNode *metaDataMT;
    float corrThresh;
} CullParams;
void cullSmooth(CullParams *cullPar);
void writeCullData(CullParams *cullPar);
void cullStats(CullParams *cullPar);
void loadCullData(CullParams *cullPar);
void loadSimData(CullParams *cullPar);
void cullSTData(CullParams *cullPar);
void cullIslands(CullParams *cullPar);
int32_t loadCullMask(CullParams *cullPar);
