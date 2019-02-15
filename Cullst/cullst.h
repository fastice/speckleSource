
typedef struct {
    int rStart;   /* Starting location for first estimate */
    int aStart;
    int deltaR;  /* Jump deltaR/A pixels between estimates */
    int deltaA;
    int nR;  /* Number of estimates in range/azimuth direction (nR*nA tot) */
    int nA;
    int maskFlag;
    int useSim;
    int singleMT;
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
    char *outFileR;    /* Range ouput file */
    char *outFileA;    /* Azimuth ouput file */
    char *outFileC;    /* Correlation ouput file */
    char *outFileT;    /* Match type ouput file */
    char *outFileD;    /* Header (.dat) file */
    char *outFileSR;    /* standard deviations */
    char *outFileSA;    
    char *inFileR;    /* Range input file */
    char *inFileA;    /* Azimuth input file */
    char *inFileC;    /* Correlation input file */
    char *inFileT;    /* Match type input file */
    char *inFileD;    /* Header (.dat) file */
    char *inFileSimR;    /* Range input file */
    char *inFileSimA;    /* Azimuth input file */
    char *inFileSimD;    /* Header (.dat) file */    

    int bR;    /* Size of window used for checking */
    int bA;
    int nGood;
    float maxA, maxR;
    int sR;    /* smoothing window size */
    int sA;
    int islandThresh;
    int ignoreOffsets;
} CullParams;
void cullSmooth(CullParams *cullPar);
void writeCullData(CullParams *cullPar);
void cullStats(CullParams *cullPar);
void loadCullData(CullParams *cullPar);
void loadSimData(CullParams *cullPar);
void cullSTData(CullParams *cullPar);
void cullIslands(CullParams *cullPar);
int loadCullMask(CullParams *cullPar);
