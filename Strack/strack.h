#include "mosaicSource/common/common.h"
#include "fft/fftw-2.1.5/fftw/fftw.h" 

#define GEODATMASK 0
#define OFFSETMASK 1

typedef short int strackComplexElement;

typedef struct strackComplexType {
	strackComplexElement r;
	strackComplexElement i;
} strackComplex;


typedef struct {
	int nr;
	int na;
	int nrl;
	int nal;
	int patchSize;
	fftwnd_plan forward;
	fftwnd_plan backward;
	fftw_complex **intf;
} intData;  

typedef struct {
	int r0;
	int a0;
	int nr;
	int na;
	int nrl;
	int nal;
	char **mask;
} maskData;  

typedef struct {
	int nr;
	int na;
	int firstRow;
	int lastRow;
	void **buf;
} StrackBuf;

               
typedef struct {
	char *imageFile1;  /* SLC image files */
	char *imageFile2;
	char *parFile1;    /* SLC par files */
	char *parFile2;    
	char *intFile;
	char *intGeodat;    /* Geodat file for inteferogram */
	char *initOffsetsFile;    /* root name for initial offsets */
	char *baseParams;  /* Baseline paramter file */
	char *baseFile;    /* Estimated Baseline file */
	char *outFileR;    /* Range ouput file */
	char *outFileA;    /* Azimuth ouput file */
	char *outFileC;    /* Correlation ouput file */
	char *outFileT;    /* Match type ouput file */
	char *outFileD;    /* Header (.dat) file */
	char *maskFile;
	char *maskGeodat;
	int maskFlag;
	int maskType;
	FILE *fpI1;        /* File pointers to images */
	FILE *fpI2;
	SARData imageP1;  /* SLC parameters */
	SARData imageP2;
	int noComplex;     /* Amp only flag */
	intData intDat;
	maskData maskDat;
	int rStart;   /* Starting location for first estimate */
	int aStart;

	int deltaR;  /* Jump deltaR/A pixels between estimates */
	int deltaA;

	int nR;  /* Number of estimates in range/azimuth direction (nR*nA tot) */
	int nA; 

	int wR;  /* Complex Window size in pixels */
	int wA;

	int wRa;  /* Amplitude Window size in pixels */
	int wAa;

	int edgePad; /* Zero pad second image by this amount */
	int edgePadA; /* Zero pad second image by this amount */
	int edgePadR; /* Zero pad second image by this amount */
	int scaleFactor; /* factor if using subsampled slc's (scales to single look)*/
	int navgR; /* Number of estimates to average internally in R/A directions*/
	int navgA; /* Number of estimates to average internally in R/A directions*/
	int intFlag; /* Flag to indicate interferogram phase will be used */
	int baseFlag; /* Flag to indicate baseline  phase will be used */
	int floatFlag; /* Flag for floating point complex input */
	int polyShift;
	/*
	  Avg power for patches
	*/
	double p1;
	double p2;
	/* Baseline Parameters */
	double Bn;
	double Bp;
	double dBn;
	double dBp;
	double dBnQ;
	double dBpQ;
	int azShift;  /* Azimuth spectral shift */
	int rangeShift; /* Range spectral shift */
	/* Offset parameters */
	char *initialOffsetFile;      /* File with shift polynomial */
	double rShiftPoly[3];  /* Rshift = [0] + [1]*rangePix + [2]*azimuthPix */
	double aShiftPoly[3];  /* Ashift = [0] + [1]*rangePix + [2]*azimuthPix */
	/* FFT plans */
	fftwnd_plan cForward;
	fftwnd_plan cForwardFast;
	fftwnd_plan cReverseNoPad;
	fftwnd_plan cReverseFast;
	fftw_plan onedForward1;
	fftw_plan onedForward2;
	fftw_plan onedForward1R;
	fftw_plan onedForward2R;
	/* Results */
	float **offR;
	float **offA;
	float **corr;
	char **type;
	/* Stats */
	int nComplex;
	int nAmp; 
	int nAmpL; 
	int nFail;
	Offsets initOffsets;
	double lambda;
	int hanningFlag;
	double latc;
} TrackParams;

void parseTrack(char *parFile, TrackParams *trackPar);
void parseBase(TrackParams *trackPar);
void parseInitialOffsets(TrackParams *trackPar);
void speckleTrack(TrackParams *trackPar);
void sTrackOut(TrackParams *trackPar);
void getInt(TrackParams *trackPar);
void getMask(TrackParams *trackPar);
