#include "mosaicSource/common/common.h"
#include "fft/fftw-2.1.5/fftw/fftw.h"

#define GEODATMASK 0
#define OFFSETMASK 1

typedef int16_t strackComplexElement;

typedef struct strackComplexType
{
	strackComplexElement r;
	strackComplexElement i;
} strackComplex;

typedef struct
{
	float x;
	float y;
} xyData;

typedef struct
{
	int32_t nr;
	int32_t na;
	int32_t nrl;
	int32_t nal;
	int32_t patchSize;
	fftwnd_plan forward;
	fftwnd_plan backward;
	fftw_complex **intf;
} intData;

typedef struct
{
	int32_t r0;
	int32_t a0;
	int32_t nr;
	int32_t na;
	int32_t nrl;
	int32_t nal;
	char **mask;
} maskData;

typedef struct
{
	int32_t nr;
	int32_t na;
	int32_t firstRow;
	int32_t lastRow;
	void **buf;
} StrackBuf;

typedef struct
{
	char *imageFile1; /* SLC image files */
	char *imageFile2;
	char *parFile1; /* SLC par files */
	char *parFile2;
	char *intFile;
	char *intGeodat;	   /* Geodat file for inteferogram */
	char *initOffsetsFile; /* root name for initial offsets */
	char *baseParams;	   /* Baseline paramter file */
	char *baseFile;		   /* Estimated Baseline file */
	char *outFileR;		   /* Range ouput file */
	char *outFileA;		   /* Azimuth ouput file */
	char *outFileC;		   /* Correlation ouput file */
	char *outFileT;		   /* Match type ouput file */
	char *outFileD;		   /* Header (.dat) file */
	char *maskFile;
	char *maskGeodat;
	int32_t maskFlag;
	int32_t maskType;
	FILE *fpI1; /* File pointers to images */
	FILE *fpI2;
	SARData imageP1; /* SLC parameters */
	SARData imageP2;
	int32_t noComplex; /* Amp only flag */
	intData intDat;
	maskData maskDat;
	int32_t rStart; /* Starting location for first estimate */
	int32_t aStart;

	int32_t deltaR; /* Jump deltaR/A pixels between estimates */
	int32_t deltaA;

	int32_t nR; /* Number of estimates in range/azimuth direction (nR*nA tot) */
	int32_t nA;

	int32_t wR; /* Complex Window size in pixels */
	int32_t wA;

	int32_t wRa; /* Amplitude Window size in pixels */
	int32_t wAa;

	int32_t edgePad;	 /* Zero pad second image by this amount */
	int32_t edgePadA;	 /* Zero pad second image by this amount */
	int32_t edgePadR;	 /* Zero pad second image by this amount */
	int32_t scaleFactor; /* factor if using subsampled slc's (scales to single look)*/
	int32_t navgR;		 /* Number of estimates to average internally in R/A directions*/
	int32_t navgA;		 /* Number of estimates to average internally in R/A directions*/
	int32_t intFlag;	 /* Flag to indicate interferogram phase will be used */
	int32_t baseFlag;	 /* Flag to indicate baseline  phase will be used */
	int32_t floatFlag;	 /* Flag for floating point32_t complex input */
	int32_t polyShift;
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
	int32_t azShift;	/* Azimuth spectral shift */
	int32_t rangeShift; /* Range spectral shift */
	/* Offset parameters */
	char *initialOffsetFile; /* File with shift polynomial */
	double rShiftPoly[3];	 /* Rshift = [0] + [1]*rangePix + [2]*azimuthPix */
	double aShiftPoly[3];	 /* Ashift = [0] + [1]*rangePix + [2]*azimuthPix */
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
	int32_t nComplex;
	int32_t nAmp;
	int32_t nAmpL;
	int32_t nFail;
	Offsets initOffsets;
	double lambda;
	int32_t hanningFlag;
	int32_t legacyFlag;
	int32_t gaussFlag;
	double latc;
	int32_t osF; /* azimuth oversample */
	int32_t maxTries;
} TrackParams;

void parseTrack(char *parFile, TrackParams *trackPar);
void parseBase(TrackParams *trackPar);
void parseInitialOffsets(TrackParams *trackPar);
void speckleTrack(TrackParams *trackPar);
void sTrackOut(TrackParams *trackPar);
void getInt(TrackParams *trackPar);
void getMask(TrackParams *trackPar);
