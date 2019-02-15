
typedef struct {
    float **XS; /* SMoothed data */ 
    float **YS;
    float **offXC; /* SMoothed data */ 
    float **offYC;
    matchResult matches;
     lsFit fitDat;
    char **type;
    char *outFileX;    /* Range ouput file */
    char *outFileY;    /* Azimuth ouput file */
    char *outFileRho;    /* Correlation ouput file */
    char *outFileT;    /* Match type ouput file */
    char *outFileD;    /* Header (.dat) file */
    char *outFileSX;    /* standard deviations */
    char *outFileSY;    
    int32 bR;    /* Size of window used for checking */
    int32 bA;
    int nGood;
    float maxX, maxY;
    int sX;    /* smoothing window size */
    int sY;
    int islandThresh;
    int32 nAttempt;
    int32 nAfterCull;
    int32 nMatch;
} CullLSParams;



void loadLSCullData(CullLSParams *cullPar);
void cullLSData(CullLSParams *cullPar);
void cullLSIslands(CullLSParams *cullPar);
void cullLSStats(CullLSParams *cullPar);
void cullLSSmooth(CullLSParams *cullPar);
