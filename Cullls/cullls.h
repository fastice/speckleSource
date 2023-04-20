
#include <stdint.h>
typedef struct
{
    float **XS; /* SMoothed data */
    float **YS;
    float **offXC; /* SMoothed data */
    float **offYC;
    matchResult matches;
    lsFit fitDat;
    char **type;
    char *outFileX;   /* Range ouput file */
    char *outFileY;   /* Azimuth ouput file */
    char *outFileRho; /* Correlation ouput file */
    char *outFileT;   /* Match type ouput file */
    char *outFileD;   /* Header (.dat) file */
    char *outFileSX;  /* standard deviations */
    char *outFileSY;
    int32_t bR; /* Size of window used for checking */
    int32_t bA;
    int32_t nGood;
    float maxX, maxY;
    int32_t sX; /* smoothing window size */
    int32_t sY;
    int32_t islandThresh;
    int32_t nAttempt;
    int32_t nAfterCull;
    int32_t nMatch;
} CullLSParams;

void loadLSCullData(CullLSParams *cullPar);
void cullLSData(CullLSParams *cullPar);
void cullLSIslands(CullLSParams *cullPar);
void cullLSStats(CullLSParams *cullPar);
void cullLSSmooth(CullLSParams *cullPar);
