#include <stdlib.h>
#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strackt.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
#include <sys/types.h>
#include "time.h"

#define NOVER 16  // Oversample for 2x oversample, (ie. 1/(NOVER*2) quantization)
#define NFAST 20  // Area to oversample for locating correlation peak
#define BAD 0
#define AMPMATCH 2
#define NBUFFERLINES 2000
#define OS 2 /* Amplitude oversample factor */
#define FFTWPLANMODE FFTW_ESTIMATE
static void correlateFast(TrackParams *trackPar, double **tmpS1, double **tmpS2, float *rShift, float *aShift, float *cMax);
static int32_t corrMatch(TrackParams *trackPar, double **tmpS1, double **tmpS2, int32_t i, int32_t j, float rShift, float aShift, float *cMax);
static int32_t getCorrPatchesFast(int32_t r1, int32_t a1, int32_t r2, int32_t a2, FILE *fp1, FILE *fp2, TrackParams *trackPar, int32_t large);
static void findImage2Pos(int32_t r1, int32_t a1, TrackParams *trackPar, int32_t *r2, int32_t *a2);
static void fftCorrPlans(TrackParams *trackPar);
static void mallocSpace(TrackParams *trackPar);
static char **mallocByteMat(int32_t nA, int32_t nR);
static float **mallocFloatMat(int32_t nx, int32_t ny);
static double **mallocDoubleMat(int32_t nA, int32_t nR);
static fftw_complex **mallocfftw_complexMat(int32_t nx, int32_t ny);
static strackComplex **mallocErs1ComplexMat(int32_t nA, int32_t nR);
static void overSampleC(TrackParams *trackPar);
static void cmpTrackFast(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax);
static void ampMatchEdge(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, int32_t large);
static void findCPeak(float **corrF, fftw_complex **corrNoPad, int32_t wR, int32_t wA, int32_t *iMax, int32_t *jMax, double *cMax);
static void writeDatFile(TrackParams *trackPar);
static void getShifts(double range, double azimuth, Offsets *offsets, double *dr, double *da);
static int32_t maskValue(TrackParams *trackPar, int32_t r1, int32_t a1);
static void initMatrix(float **data,  int32_t wA, int32_t wR, float value);
static void initByteMatrix(char **data,  int32_t wA, int32_t wR, char value);
static void zeroPad(fftw_complex **f1, fftw_complex **f2, fftw_complex **f1a, fftw_complex **f2a, TrackParams *trackPar);
static void getComplexData(int32_t a1, int32_t a2, int32_t r1, int32_t r2, fftw_complex **im1in, fftw_complex **im2in, TrackParams *trackPar);
static void detectPatch(float **data, fftw_complex **im, int32_t wR, int32_t wA, int32_t edgePadR, 
						int32_t edgePadA, float scale, TrackParams *trackPar, int32_t largePatch);
static void updateSLCBuffer(int32_t bufferNum, TrackParams *trackPar, int32_t a1, int32_t sSize, FILE *fp);
static void getPeakCorr(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double meanR, double sigmaR, double *maxCorr);
static void getInitialGuess(TrackParams *trackPar);
/* ************************ GLOBAL VARIABLES ***************************/
fftw_complex **patch1in, **patch2in; /* ... */
fftw_complex **psNoPad;
fftw_complex **cFast, **cFastOver;
fftw_complex **cNoPad;
fftw_complex **psFast, **psFastOver;
time_t startTime, lastTime;
/* Amplitude stuff */
fftw_complex **img1, **img2; /* Detected images for amplitude match*/
fftw_complex **psAmpNoPad;
fftw_complex **fftFa1, **fftFa2; /* FFt's */
fftw_complex **caNoPad;
float **caNoPadMag;
/* FFT plans */
fftwnd_plan aForward;
fftwnd_plan aReverseNoPad;
StrackBuf imageBuf1;
StrackBuf imageBuf2;
/* New stuff for oversampling ^^^ */
fftw_complex **img1in, **img2in; /* ^^^ Detected images for amplitude match*/
fftwnd_plan aForwardIn;
fftw_complex **fftFa1os, **fftFa2os; 
fftw_complex **patch1, **patch2;
double **meanS;
double **sigmaS;
double **corrResult;
float **dataS;
float **dataR;

/*
  Main routine for amplitude matching with an edge pad.
*/

void corrTrackFast(TrackParams *trackPar)
{
	//extern float **cFastOverMag; //, **cNoPadMag;
	extern fftw_complex **psNoPad;
	extern double **meanS;
	extern double **sigmaS;
	extern double **corrResult;
	double **tmpS1, **tmpS2;
	extern float **dataS, **dataR;
	FILE *fp1, *fp2;
	FILE *fpR, *fpA, *fpC, *fpT;
	int32_t i, j, nGood;
	int32_t r1, a1, r2, a2;
	int32_t wR2, wA2;
	float cMax;
	int32_t WA2, WR2;
	extern time_t startTime, lastTime;
	float rShift, aShift;
	int32_t nTot;
	int32_t good;
	double cAvg;
	int32_t maskVal; /* Flag to determine whether to match */
	int32_t nMask;
	fprintf(stderr, "\nSPECKLE TRACKING\n");
	/*
	   Internally keep edge pad bigger to acount for over sampling around the peak.
	   Program operates though edge pad is size specified (e.g., peak must be within
	   original edge pad limits.
	*/
	trackPar->edgePadR += NFAST / 4;
	trackPar->edgePadA += NFAST / 4;
	// Malloc space
	tmpS1 = mallocDoubleMat(trackPar->wAa * OS, trackPar->wRa * OS);
	tmpS2 = mallocDoubleMat(trackPar->wAa * OS, trackPar->wRa * OS);
	mallocSpace(trackPar);
	// Get size of smaller search window
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	fprintf(stderr, "FFT Size (not oversampled) R,A %i %i\n", trackPar->wRa, trackPar->wAa);
	fprintf(stderr, "Search chip  (not oversampled) R,A %i %i\n", trackPar->wRa - 2 * trackPar->edgePadR, trackPar->wAa - 2 * trackPar->edgePadA);
	fprintf(stderr, "Search Radius %i %i\n", trackPar->edgePadR - NFAST / 4, trackPar->edgePadA - NFAST / 4);
	fprintf(stderr, "Search Radius with PAD %i %i\n", trackPar->edgePadR, trackPar->edgePadA);
	fprintf(stderr, "wR2,wA2 %i %i\n", wR2 / 2, wA2 / 2);
	// Open output files
	fpR = fopen(trackPar->outFileR, "w");
	fpA = fopen(trackPar->outFileA, "w");
	fpC = fopen(trackPar->outFileC, "w");
	fpT = fopen(trackPar->outFileT, "w");
	// Open input files
	fp1 = openInputFile(trackPar->imageFile1);
	fp2 = openInputFile(trackPar->imageFile2);
	// FFTW Plan
	fprintf(stderr, "Plans\n");
	fftCorrPlans(trackPar);
	// Initializations
	trackPar->nFail = 0.0;
	trackPar->nComplex = 0.0;
	trackPar->nAmp = 0;
	nMask = 0;
	// Initial offsets initialization
	getInitialGuess(trackPar);
	// Loop on azimuth
	for (i = 0; i < trackPar->nA; i++)
	{
		lastTime = time(NULL);
		// Azimuth coordinate
		a1 = trackPar->aStart + i * trackPar->deltaA - trackPar->wAa / 2;
		cAvg = 0; // Avg corr for the row
		nGood = 0; // nGood for the row
		// Loop on range 
		for (j = 0; j < trackPar->nR; j++)
		{
			r1 = trackPar->rStart + j * trackPar->deltaR - trackPar->wRa / 2;
			// Find position for second image
			findImage2Pos(r1, a1, trackPar, &r2, &a2);
			rShift = r2 - r1;
			aShift = a2 - a1;
			cMax = 0.;
			/* Fixed 10/29/18 - corrections to r1/a1 applied in maskValue */
			maskVal = maskValue(trackPar, r1 * trackPar->scaleFactor, a1 * trackPar->scaleFactor);
			if (maskVal != 0 && maskVal != 1)
				fprintf(stderr, "maskVal %i\n", maskVal);
			nMask += maskVal;
			// If indicated my mask value, do the match
			if (maskVal == 1)
			{
				/* check bounds of start postion for each patch (lower left corner) */
				if (r1 > 0 && r1 < trackPar->imageP1.nSlpR && a1 > 0 && a1 < trackPar->imageP1.nSlpA &&
					r2 > 0 && r2 < trackPar->imageP2.nSlpR && a2 > 0 && a2 < trackPar->imageP2.nSlpA)
				{
					/* Read patches for amplitude match */
					getCorrPatchesFast(r1, a1, r2, a2, fp1, fp2, trackPar, FALSE);
					/* Do the  match, load values to trackPar->offR/A/corr/type */
					corrMatch(trackPar, tmpS1, tmpS2, i, j, rShift, aShift, &cMax);
				}
				else
					cMax = 0.0; /* end if r1>0...*/
				// Update stats
				if (cMax > .00)
				{
					cAvg += cMax;
					nGood++;
					good = TRUE;
				} else good = FALSE;

			}
			else
			{
				trackPar->offR[i][j] = -LARGEINT;
				trackPar->offA[i][j] = -LARGEINT;
				trackPar->corr[i][j] = 0;
				trackPar->type[i][j] = 0;
			} /* end if maskVal..*/
		} // End range loop
		writeOffsets(i, trackPar, fpR, fpA, fpC, fpT);
		// Update stats
		nTot = (i + 1) * trackPar->nR;
		if (nGood > 0)
			cAvg = cAvg / (double)(nGood);
		else
			cAvg = 0;
		fprintf(stderr,
				"\r%6i nTot %i, nMatch %8i %4.1f  "
				"nAmp %8i %4.1f (%7i) nFail %7i %4.1f cAvg(line) %4.2f %i -- %5i (s) --",
				a1, nTot, trackPar->nAmp, 100. * (double)(trackPar->nAmp ) / (double)nTot,
				trackPar->nAmp ,
				100. * (double)(trackPar->nAmp) / (double)(nTot),
				trackPar->nAmp, trackPar->nFail,
				100. * (double)(trackPar->nFail) / (double)(nTot),
				cAvg, nMask, (int)(time(NULL) - lastTime));
	} /* End for i=0... */
	fclose(fpR);
	fclose(fpA);
	fclose(fpC);
	fclose(fpT);
	writeDatFile(trackPar);
	writeVrtFile(trackPar);
}

// Find inital shift
  
static void findImage2Pos(int32_t r1, int32_t a1, TrackParams *trackPar, int32_t *r2, int32_t *a2)
{
	double aShift, rShift;
	double sc;
	sc = (double)trackPar->scaleFactor;
	if (trackPar->polyShift == TRUE)
	{
		rShift = trackPar->rShiftPoly[0] + (double)(r1 * sc) * trackPar->rShiftPoly[1] +
				 (double)(a1 * sc) * trackPar->rShiftPoly[2];
		aShift = trackPar->aShiftPoly[0] + (double)(r1 * sc) * trackPar->aShiftPoly[1] +
				 (double)(a1 * sc) * trackPar->aShiftPoly[2];
	}
	else
	{
		getShifts(r1, a1, &(trackPar->initOffsets), &rShift, &aShift);
		/* If offset not available use poly - added 6/2/04 */
		if (rShift < (-LARGEINT + 1))
		{
			rShift = trackPar->rShiftPoly[0] + (double)(r1 * sc) * trackPar->rShiftPoly[1] +
					 (double)(a1 * sc) * trackPar->rShiftPoly[2];
			aShift = trackPar->aShiftPoly[0] + (double)(r1 * sc) * trackPar->aShiftPoly[1] +
					 (double)(a1 * sc) * trackPar->aShiftPoly[2];
		}
	}
	/* Added 8/21 for tsx scaling */
	rShift /= sc;
	aShift /= sc;
	if (rShift >= 0)
		*r2 = r1 + (int)(rShift + 0.5);
	else
		*r2 = r1 + (int)(rShift - 0.5);
	if (aShift >= 0)
		*a2 = a1 + (int)(aShift + 0.5);
	else
		*a2 = a1 + (int)(aShift - 0.5);
}

//  Read patches
//  It reads the full size patch for image 2, so that subsequent programs
//  will just read the relevant patch when needed.
static int32_t getCorrPatchesFast(int32_t r1, int32_t a1, int32_t r2, int32_t a2,
								  FILE *fp1, FILE *fp2, TrackParams *trackPar, int32_t large)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	extern fftw_complex **img1, **img2;
	extern fftw_complex **img1in, **img2in;	// Detected images for amplitude match
	extern fftwnd_plan aForwardIn;
	extern fftw_complex **fftFa1os, **fftFa2os;
	extern fftw_complex **fftFa1, **fftFa2; /* FFt's */
	int32_t wR2, wA2, wRa, wAa;
	size_t sSize;
	float scale;
	if (a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA)
	{
		error("getCorrPatchesFast");
		return (FALSE);
	}
	/* set points to buffer for image 1 and image 2*/
	wAa = trackPar->wAa;
	wRa = trackPar->wRa;
	if (trackPar->floatFlag == TRUE)
		sSize = sizeof(fftw_complex);
	else
		sSize = sizeof(strackComplex);
	/* size for search chip */
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	// Added navg values to scale computation to avoid overflow that was 
	//  occurring with tsx data 
	scale = 1. / ((float)wRa * (float)wAa * (float)wRa * (float)wAa *
				  (float)(OS * OS * OS * OS) * (OS * OS * (trackPar->navgA + 1) * (trackPar->navgR + 1)));		  
	//  Load buffers if needed (only required if patch is outside buffer)
	updateSLCBuffer(1, trackPar, a1, sSize, fp1);
	updateSLCBuffer(2, trackPar, a2, sSize, fp2);
	// Get complex data
	getComplexData(a1, a2, r1, r2, img1in, img2in, trackPar);
	// Forward FFT 
	fftwnd_one(aForwardIn, img1in[0], fftFa1os[0]);
	fftwnd_one(aForwardIn, img2in[0], fftFa2os[0]);
	// Zero pad fft for over sampling
	zeroPad(fftFa1os, fftFa2os, fftFa1, fftFa2, trackPar);
	// Inverse transform
	fftwnd_one(aReverseNoPad, fftFa1[0], img1[0]);
	fftwnd_one(aReverseNoPad, fftFa2[0], img2[0]);
	//  Detect Data , Note the image patches are double bookkept as im1/im2 and dataR, dataS
	//  based on the way the program was kluged together
	// Smaller patch
	detectPatch(dataR, img2, wR2, wA2, trackPar->edgePadR, trackPar->edgePadA, scale, trackPar, FALSE);
	// Larger patch
	detectPatch(dataS, img1, trackPar->wRa * OS, trackPar->wAa * OS, 0, 0, scale, trackPar, TRUE);
	/*
	FILE *fpDebug;
	fprintf(stderr, "%li\n", trackPar->wRa *OS * trackPar->wAa *OS*sizeof(fftw_complex));
	fpDebug= fopen("test1.debug","w");	fwriteBS(dataS[0],trackPar->wRa *OS * trackPar->wAa *OS,sizeof(float),fpDebug,FLOAT32FLAG); fclose(fpDebug);
	fpDebug= fopen("test2.debug","w");	fwriteBS(dataR[0],wR2*wA2, sizeof(float), fpDebug,FLOAT32FLAG); fclose(fpDebug);
	fpDebug= fopen("test1a.debug","w");	fwriteBS(img1[0],trackPar->wRa *OS * trackPar->wAa *OS ,sizeof(fftw_complex), fpDebug,FLOAT32FLAG); fclose(fpDebug);
	fpDebug= fopen("test2a.debug","w");	fwriteBS(img2[0],trackPar->wRa *OS * trackPar->wAa *OS,sizeof(fftw_complex),fpDebug,FLOAT32FLAG); fclose(fpDebug);
		 exit(-1);
	*/
	return (TRUE);
}


/***************************************************************
	Main driver for amplitude match (with edgepad)
***************************************************************/
static int32_t corrMatch(TrackParams *trackPar, double **tmpS1, double **tmpS2, int32_t i, int32_t j, float rShift, float aShift, float *cMax)
{
	extern float **dataS, **dataR;
	float rShift1, aShift1;
	int32_t wR2, wA2;
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	// Do the correlation
	correlateFast(trackPar, tmpS1, tmpS2, &rShift1, &aShift1, cMax);
	// Accept or reject match
	if (*cMax > .00)
	{
		rShift -= rShift1; // Remove the initial guess
		aShift -= aShift1;
		trackPar->nAmp++;
	}
	else
	{
		rShift = -LARGEINT;
		aShift = -LARGEINT;
		trackPar->nFail++;
	}
	// Save values 
	trackPar->offR[i][j] = rShift;
	trackPar->offA[i][j] = aShift;
	trackPar->corr[i][j] = *cMax;
	trackPar->type[i][j] = AMPMATCH; 
	return TRUE;
}

// Compute mean and sigma. Return values for small patch, large patch is saved for each position in meanS, sigmaS
static void computeMeanSigma(TrackParams *trackPar, double *meanR,  double *sigmaR, double **tmpS1, double **tmpS2) 
{
	extern float **dataS, **dataR;
	extern double **meanS;
	extern double **sigmaS;
	double sum1, sum2;
	int32_t lr, la, s, l, s1, l1, s2, l2;
	int32_t wR2, wA2;
	*sigmaR = 0;
	*meanR = 0;
	// Compute single values for small chip
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	for (la = 0; la < wA2; la++) 
	{
		for (lr = 0; lr < wR2; lr++)
		{
			*sigmaR += dataR[la][lr] * dataR[la][lr];
			*meanR += dataR[la][lr];
		}
	}
	*meanR /= (double)(wR2 * wA2);
	*sigmaR = *sigmaR / (double)(wR2 * wA2 - 1.);
	*sigmaR -= *meanR * (*meanR);
	//   Step 2:
	//   Compute mean and variance for large image for each 
	//  potential chip position. First,do running average over columns
	for (l = 0; l < trackPar->wAa * OS; l++)
	{
		sum1 = 0.0;
		sum2 = 0.0;
		s1 = 0;
		for (s = 0; s < trackPar->wRa * OS; s++)
		{
			sum1 += dataS[l][s];
			sum2 += dataS[l][s] * dataS[l][s];
			if (s >= (wR2 - 1))
			{
				s2 = s - (wR2 - 1);
				tmpS1[l][s1] = sum1;
				tmpS2[l][s1] = sum2;
				sum1 -= dataS[l][s2];
				sum2 -= dataS[l][s2] * dataS[l][s2];
				s1++;
			}
		}
	}
	// finish by doing running average over rows
	s1 = 0;
	for (s = 0; s < (2 * trackPar->edgePadR * OS) + 1; s++)
	{
		sum1 = 0.0;
		sum2 = 0.0;
		l1 = 0;
		for (l = 0; l < trackPar->wAa * OS; l++)
		{
			sum1 += tmpS1[l][s];
			sum2 += tmpS2[l][s];
			if (l >= (wA2 - 1))
			{
				l2 = l - (wA2 - 1);
				meanS[l1][s1] = sum1 / (double)(wR2 * wA2);
				sigmaS[l1][s1] = sum2 / (double)(wR2 * wA2 - 1) - meanS[l1][s1] * meanS[l1][s1];
				sum1 -= tmpS1[l2][s];
				sum2 -= tmpS2[l2][s];
				l1++;
			}
		}
		s1++;
	}
}

/***********************************************************
   Main code to execute correlation
************************************************************/

static void correlateFast(TrackParams *trackPar, double **tmpS1, double **tmpS2,  float *rShift, float *aShift, float *cMax)
{
	extern float **dataS, **dataR;
	extern double **meanS, **sigmaS;
	extern fftw_complex **cFast, **cFastOver, **caNoPad;
	double maxCorr, corr;
	extern double **corrResult;
	int32_t iMax, jMax, iMax1, jMax1;
	double rShift1, aShift1;
	double cM;
	int32_t i, j, i1, j1;
	int32_t wR2, wA2;
	double meanR, sigmaR;
	
	double scaleW;
	// size of search chip 
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	// Step 1: Compute mean/sigma for second image patch
	computeMeanSigma(trackPar, &meanR,  &sigmaR, tmpS1, tmpS2);
	// Step 3, do fft convolution of image patchtes
	ampMatchEdge(trackPar, &iMax, &jMax, &cM, FALSE);
	// Step 4: find peak in correlation function
	getPeakCorr(trackPar, &iMax, &jMax, meanR, sigmaR, &maxCorr);
	// Return if no data
	if (maxCorr < 0)
		return;
	// Step 5:  load from around peak to oversample
	for (i = 0; i <= NFAST; i++)
	{
		i1 = iMax - NFAST / 2 + i;
		j1 = jMax - NFAST / 2;
		for (j = 0; j <= NFAST; j++)
		{
			cFast[i][j].re = corrResult[i1][j1];
			cFast[i][j].im = 0.0;
			j1++;
		}
	}
	// Step 6: Oversample correlation peak
	overSampleC(trackPar);
	// Step 7: Find peak and renorm value
	*cMax = -1.0e30;
	for (i = ((NFAST+1) * NOVER) / 2 - NOVER * 2; i < ((NFAST+1) * NOVER) / 2 + NOVER * 2; i++)
		for (j = ((NFAST+1) * NOVER) / 2 - NOVER * 2; j < ((NFAST+1) * NOVER) / 2 + NOVER * 2; j++)
		{
			corr = (cFastOver[i][j].re * cFastOver[i][j].re + cFastOver[i][j].im * cFastOver[i][j].im);
			if (corr > *cMax)
			{
				*cMax = corr;
				iMax1 = i;
				jMax1 = j;
			}
		}
	*cMax = sqrt(*cMax) / ( (NFAST+1) * (NFAST+1));	
	// Step 8: Compute raw, fractional pixel shift.
	iMax = (iMax * NOVER) + iMax1 - (NOVER * (NFAST)) / 2;
	jMax = (jMax * NOVER) + jMax1 - (NOVER * (NFAST)) / 2;
	*rShift = ((float)jMax - trackPar->edgePadR * OS * NOVER) / (NOVER * OS);
	*aShift = ((float)iMax - trackPar->edgePadA * OS * NOVER) / (NOVER * OS);
}

static void updateSLCBuffer(int32_t bufferNum, TrackParams *trackPar, int32_t a1, int32_t sSize, FILE *fp)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	StrackBuf *imageBuf;
	int32_t nSlpA, nSlpR, a1a, s1;
	off_t offset1;
	// Sellect buffer
	if (bufferNum == 1)
	{
		imageBuf = &imageBuf1;
		nSlpA = trackPar->imageP1.nSlpA;
		nSlpR = trackPar->imageP1.nSlpR;
	}
	else if (bufferNum == 2)
	{
		imageBuf = &imageBuf2;
		nSlpA = trackPar->imageP2.nSlpA;
		nSlpR = trackPar->imageP2.nSlpR;
	}
	else
		error("INVALID BUFFER");
	// load buffer
	if (a1 < imageBuf->firstRow || (a1 + trackPar->wAa) > imageBuf->lastRow)
	{
		a1a = max(0, a1 - trackPar->wAa);
		/*a1a=max(0,a1-trackPar->wRa);****/
		s1 = nSlpR * min(nSlpA - a1a, NBUFFERLINES);
		offset1 = ((off_t)a1a * nSlpR) * (off_t)sSize;
		fseeko(fp, offset1, SEEK_SET);
		if (trackPar->floatFlag == TRUE)
			freadBS((void *)imageBuf->buf[0], sSize, s1, fp, FLOAT32FLAG);
		else
			freadBS((void *)imageBuf->buf[0], sSize, s1, fp, INT16FLAG);
		imageBuf->firstRow = a1a;
		imageBuf->lastRow = a1a + min(nSlpA - a1a, NBUFFERLINES) - 1;
	}
}

// Mv data from buffer to input image patches
static void getComplexData(int32_t a1, int32_t a2, int32_t r1, int32_t r2, fftw_complex **im1in,
						   fftw_complex **im2in, TrackParams *trackPar)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	size_t sSize;
	int32_t i1, i2, i, j, j1, j2;
	strackComplex **ers1Buf1, **ers1Buf2;
	fftw_complex **fftwBuf1, **fftwBuf2;
	// Set up for float or int data
	if (trackPar->floatFlag == TRUE)
	{
		sSize = sizeof(fftw_complex);
		fftwBuf1 = (fftw_complex **)imageBuf1.buf;
		fftwBuf2 = (fftw_complex **)imageBuf2.buf;
	}
	else
	{
		sSize = sizeof(strackComplex);
		ers1Buf1 = (strackComplex **)imageBuf1.buf;
		ers1Buf2 = (strackComplex **)imageBuf2.buf;
	}
	// Load data from image buffers
	for (i = 0; i < trackPar->wAa; i++)
	{
		i1 = a1 - imageBuf1.firstRow + i;
		i2 = a2 - imageBuf2.firstRow + i;
		if (i1 > NBUFFERLINES || i2 > NBUFFERLINES)
			error("Buffer exceeded in getComplex Data \n");
		for (j = 0; j < trackPar->wRa; j++)
		{
			j1 = r1 + j;
			j2 = r2 + j;
			if (trackPar->floatFlag == TRUE)
			{
				im1in[i][j].re = (fftw_real)fftwBuf1[i1][j1].re;
				im1in[i][j].im = (fftw_real)fftwBuf1[i1][j1].im;
				im2in[i][j].re = (fftw_real)fftwBuf2[i2][j2].re;
				im2in[i][j].im = (fftw_real)fftwBuf2[i2][j2].im;
			}
			else
			{
				im1in[i][j].re = ((fftw_real)ers1Buf1[i1][j1].r);
				im2in[i][j].re = ((fftw_real)ers1Buf2[i2][j2].r);
				im1in[i][j].im = ((fftw_real)ers1Buf1[i1][j1].i);
				im2in[i][j].im = ((fftw_real)ers1Buf2[i2][j2].i);
			}
		}
	}
}

// Zero pad for oversampling patches
static void zeroPad(fftw_complex **f1, fftw_complex **f2, fftw_complex **f1a, fftw_complex **f2a, TrackParams *trackPar)
{
	int32_t i, i1, i2, ia, j, j1, j2, j2a, ja, wRa, wAa;
	fftw_complex zero;
	zero.re = 0;
	zero.im = 0;
	wRa = trackPar->wRa;
	wAa = trackPar->wAa;
	// Zero the whole ff
	for (i = 0; i < wAa * OS; i++)
	{
		for (j = 0; j < wRa * OS; j++)
		{
			f1a[i][j] = zero;
			f2a[i][j] = zero;
		}
	}
	// Fill the corners
	for (i = 0; i < wAa / 2; i++)
	{
		i1 = OS * wAa - wAa / OS + i;
		i2 = (wAa / OS + i) % wAa;
		ia = i % wAa;
		j1 = wRa * OS - wRa / OS;
		j2 = wRa / OS;
		for (j = 0; j < wRa / 2; j++)
		{
			ja = j % wRa;
			j2a = j2 % wRa;
			f1a[i][j] = f1[ia][ja];
			f1a[i1][j] = f1[i2][ja];
			f1a[i1][j1] = f1[i2][j2a];
			f1a[i][j1] = f1[ia][j2a];
			f2a[i][j] = f2[ia][ja];
			f2a[i1][j] = f2[i2][ja];
			f2a[i1][j1] = f2[i2][j2a];
			f2a[i][j1] = f2[ia][j2a];
			j1++;
			j2++;
		}
	}
}

static void detectPatch(float **data, fftw_complex **im, int32_t wR, int32_t wA, int32_t edgePadR,
						int32_t edgePadA, float scale, TrackParams *trackPar, int32_t largePatch)
{
	int32_t i, j, i1, j1, k, l, n, m, extraR, extraA;
	fftw_complex zero;
	zero.re = 0.0;
	zero.im = 0.0;
	// Small catch can exceed window, but not the large patch
	if(largePatch == TRUE) 
	{
		extraR = 0;
		extraA = 0;
	}
	else
	{
		extraR = trackPar->navgR; 
		extraA = trackPar->navgA;
	}
	for (i = 0; i < wA; i++)
	{
		i1 = edgePadA * OS + i;
		for (j = 0; j < wR; j++)
		{
			j1 = edgePadR * OS + j;
			// Smoothing
			if (trackPar->navgA > 1 || trackPar->navgR > 1)
			{
				data[i][j] = 0;
				for (m = -trackPar->navgA; m <= trackPar->navgA; m++)
				{
					k = min(max(i1 + m, 0), edgePadR * OS + wA + extraA - 1); /* Changed 11/12/08 data outside window */
					for (n = -trackPar->navgR; n <= trackPar->navgR; n++)
					{
						/* fixed 11/12/08, changed from m to n */
						l = min(max(j1 + n, 0), edgePadR * OS + wR + extraR - 1);
						data[i][j] += (im[k][l].re * im[k][l].re + im[k][l].im * im[k][l].im) * scale;
					}
				}
			}
			else
			{
				data[i][j] = (im[i1][j1].re * im[i1][j1].re + im[i1][j1].im * im[i1][j1].im) * scale;
			}
		}
	}
	// Copy to im with second loop to avoid smoothing issues.
	for (i = 0; i < wA; i++)
	{
		i1 = edgePadA * OS + i;
		for (j = 0; j < wR; j++)
		{
			j1 = edgePadR * OS + j;
			data[i][j] = sqrt(data[i][j]);
			im[i1][j1].re = data[i][j];
			im[i1][j1].im = 0.0;
		}
	}
	// Zero edges if there is an edgpad
	for (i = 0; i < trackPar->wAa * OS; i++)
	{
		for (j = 0; j < edgePadR * OS; j++)
		{
			im[i][j] = zero;
			im[i][trackPar->wRa * OS - j - 1] = zero;
		}
	}
	for (i = 0; i < edgePadA * OS; i++)
	{
		for (j = 0; j < trackPar->wRa * OS; j++)
		{
			im[i][j] = zero;
			im[trackPar->wAa * OS - i - 1][j] = zero;
		}
	}
}

// Interp offsets used for initial guess
static float interpOffsets(float **image, float range, float azimuth, float minvalue, 
						   int32_t nr, int32_t na)
{
	float p1, p2, p3, p4, t, u, dr;
	int32_t i, j;
	j = (int)range;
	i = (int)azimuth;
	// out of range
	if (range < 0.0 || azimuth < 0.0 || (int)range >= nr || (int)azimuth >= na)
		return -LARGEINT;
	// Handle border pixels
	if (j == (nr - 1) || i == (int)(na - 1))
	{
		dr = image[i][j];
		if (dr < minvalue)
			return -LARGEINT;
	}
	p1 = image[i][j];
	p2 = image[i][j + 1];
	p3 = image[i + 1][j + 1];
	p4 = image[i + 1][j];
	if (p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue)
	{
		dr = max(max(max(p1, p2), p3), p4);
		if (dr < minvalue)
			return -LARGEINT;
		return dr;
	}
	return (float)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 +
				    t * u * p3 + (1.0 - t) * u * p4);
}

/***************************************************************************
   Get shifts from previous estimate
****************************************************************************/
static void getShifts(double range, double azimuth, Offsets *offsets, double *dr, double *da)
{
	float minvalue;

	minvalue = (float)-LARGEINT * 0.9;
	range = (range - offsets->rO) / offsets->deltaR;
	azimuth = (azimuth - offsets->aO) / offsets->deltaA;
	/* Out of bounds */
	*dr = interpOffsets((float **)offsets->dr, range, azimuth, minvalue, offsets->nr, offsets->na);
	*da = interpOffsets((float **)offsets->da, range, azimuth, minvalue, offsets->nr, offsets->na);
	if (*dr <= minvalue || *da <= minvalue)
	{
		*da = minvalue;
		*dr = minvalue;
	}
}

// Find the peak of the initial correlation function
static void getPeakCorr(TrackParams *trackPar, int32_t *iMax, int32_t *jMax,
					    double meanR, double sigmaR, double *maxCorr)
{
	double scaleW;
	extern double **meanS;
	extern double **sigmaS;
	int32_t i, j, l, m, n, s;
	double wR2, wA2;
	wR2 = trackPar->wRa - 2 * trackPar->edgePadR;
	wA2 = trackPar->wAa - 2 * trackPar->edgePadA;
	scaleW = 1.0 / ((OS * OS * OS * OS) * trackPar->wAa * wA2 * trackPar->wRa * wR2);
	*maxCorr = -1e30;
	// Loop to search over.Note goes past edge pad to fill buffer for over sampling */
	for (l = -OS * trackPar->edgePadA; l <= OS * trackPar->edgePadA; l++)
	{
		i = l + OS * trackPar->edgePadA;
		m = l + trackPar->wAa * OS / 2;
		for (s = -OS * trackPar->edgePadR; s <= OS * trackPar->edgePadR; s++)
		{
			j = s + OS * trackPar->edgePadR;
			n = s + trackPar->wRa * OS / 2;
			corrResult[i][j] = caNoPadMag[m][n] * scaleW;
			corrResult[i][j] = corrResult[i][j] - (meanR * meanS[i][j]);
			corrResult[i][j] = corrResult[i][j] / sqrt(sigmaR * sigmaS[i][j]);
			// Only select point32_t if within original specified edge pad 
			if (corrResult[i][j] > *maxCorr &&
				l >= -OS * (trackPar->edgePadA - NFAST / 4) &&
				s >= -OS * (trackPar->edgePadR - NFAST / 4) &&
				l <= OS * (trackPar->edgePadA - NFAST / 4) &&
				s <= OS * (trackPar->edgePadR - NFAST / 4))
			{
				*maxCorr = corrResult[i][j];
				*iMax = i;
				*jMax = j;
			}
		}
	}
}


//   FFT's to do cross correlation for matching
static void ampMatchEdge(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, int32_t large)
{
	extern fftwnd_plan aForward;
	extern fftwnd_plan aReverseNoPad;
	extern fftw_complex **psAmpNoPad;
	extern float **caNoPadMag;
	extern fftw_complex **caNoPad;
	extern fftw_complex **fftFa1, **fftFa2;
	extern fftw_complex **img1, **img2;
	int32_t wR, wA;
	int32_t i, j, i1, j1;
	wA = trackPar->wAa;
	wR = trackPar->wRa;
	// Step 1: FFT images  for cross-power spectrum
	fftwnd_one(aForward, img1[0], fftFa1[0]);
	fftwnd_one(aForward, img2[0], fftFa2[0]);
	// Step 2: Compute zero padded power spectrum
	i1 = 0;
	for (i = 0; i < wA * OS; i++)
	{
		for (j = 0; j < wR * OS; j++)
		{
			psAmpNoPad[i][j].re = fftFa1[i1][j].re * fftFa2[i1][j].re + fftFa1[i1][j].im * fftFa2[i1][j].im;
			psAmpNoPad[i][j].im = fftFa1[i1][j].im * fftFa2[i1][j].re - fftFa1[i1][j].re * fftFa2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
	// Step 3: Inverse transform to get convolved image patches
	fftwnd_one(aReverseNoPad, psAmpNoPad[0], caNoPad[0]);
	// Step 4: Unscramble result
	for (i = 0; i < wA * OS / 2; i++)
	{
		i1 = (wA * OS) / 2 + i;
		j1 = (wR * OS) / 2;
		for (j = 0; j < wR * OS / 2; j++)
		{
			// Compute magnitude
			caNoPadMag[i][j] = sqrt(caNoPad[i1][j1].re * caNoPad[i1][j1].re + caNoPad[i1][j1].im * caNoPad[i1][j1].im);
			caNoPadMag[i1][j1] = sqrt(caNoPad[i][j].re * caNoPad[i][j].re + caNoPad[i][j].im * caNoPad[i][j].im);
			caNoPadMag[i1][j] = sqrt(caNoPad[i][j1].re * caNoPad[i][j1].re + caNoPad[i][j1].im * caNoPad[i][j1].im);
			caNoPadMag[i][j1] = sqrt(caNoPad[i1][j].re * caNoPad[i1][j].re + caNoPad[i1][j].im * caNoPad[i1][j].im);
			j1++;
		}
	}
}

/*************************************************************
   Over sample correlation
**************************************************************/
static void overSampleC(TrackParams *trackPar)
{
	extern fftw_complex **psFast, **psFastOver;
	extern fftw_complex **cFast, **cFastOver;
	int32_t i, j;
	int32_t i1, i2, j1, j2;
	int32_t half;
	half = NFAST/2;
	fftwnd_one(trackPar->cForwardFast, cFast[0], psFast[0]);
	//  Zero pad fft for over sampling
	for (i = 0; i <= half; i++)
	{
		i1 = (NFAST + 1) * NOVER - half - 1 + i;
		i2 = half  + i;
		j1 = (NFAST + 1) * NOVER - half - 1;
		j2 = half;
		for (j = 0; j <= half; j++)
		{	// lower left corner
			psFastOver[i][j].re = psFast[i][j].re;
			psFastOver[i][j].im = psFast[i][j].im;
			if(j != 0) { // Upper left corner
				psFastOver[i][j1].re = psFast[i][j2].re;
				psFastOver[i][j1].im = psFast[i][j2].im;
			}	
			if(i != 0) { 
				// lower right corner
				psFastOver[i1][j].re = psFast[i2][j].re;
				psFastOver[i1][j].im = psFast[i2][j].im;
				if( j != 0 ) 
				{   // upper right corner
					psFastOver[i1][j1].re = psFast[i2][j2].re;
					psFastOver[i1][j1].im = psFast[i2][j2].im;
				}	
			}
			j1++;
			j2++;
		} /* End for j */
	} /* End for i */
	fftwnd_one(trackPar->cReverseFast, psFastOver[0], cFastOver[0]);
}

/******************************************************************************
	Do plans for fft's
	************************************************************** */

static void fftCorrPlans(TrackParams *trackPar)
{
	extern fftwnd_plan aForwardIn, aForward;
	extern fftwnd_plan aReverseNoPad;
	fprintf(stderr, "corrPatch %i %i\n", NFAST + 1, NFAST+1);
	trackPar->cForwardFast = fftw2d_create_plan(NFAST + 1, NFAST +1 , FFTW_FORWARD, FFTWPLANMODE);
	fprintf(stderr, "corrOverPatch %i %i\n", NOVER * (NFAST+1), NOVER * (NFAST+1));
	trackPar->cReverseFast = fftw2d_create_plan(NOVER * (NFAST+1), NOVER * (NFAST+1), FFTW_BACKWARD, FFTWPLANMODE);
	fprintf(stderr, "5 %i\n", trackPar->wA);
	fprintf(stderr, "6- %i %i\n", trackPar->wAa * OS, trackPar->wRa * OS);
	aForward = fftw2d_create_plan(trackPar->wAa * OS, trackPar->wRa * OS, FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "8 - %i %i \n", trackPar->wAa * OS, trackPar->wRa * OS);
	aReverseNoPad = fftw2d_create_plan(trackPar->wAa * OS, trackPar->wRa * OS, FFTW_BACKWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "10 - %i %i\n", trackPar->wAa, trackPar->wRa);
	aForwardIn = fftw2d_create_plan(trackPar->wAa, trackPar->wRa, FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
}

static void getInitialGuess(TrackParams *trackPar)
{
if (trackPar->polyShift == FALSE)
	{
		if (readBothOffsetsStrackVrt(&(trackPar->initOffsets), trackPar->initOffsetsFile) == FALSE)
		{
			// Old method
			trackPar->initOffsets.file = appendSuff(trackPar->initOffsetsFile, ".da",
													malloc(strlen(trackPar->initOffsetsFile) + 4));
			trackPar->initOffsets.rFile = appendSuff(trackPar->initOffsetsFile, ".dr",
													 malloc(strlen(trackPar->initOffsetsFile) + 4));
			readBothOffsetsStrack(&(trackPar->initOffsets));
		}
	}
}

/*************************************************************
   Malloc space for all the global arrays/matrices and elements of trackPar
*************************************************************/

static void mallocSpace(TrackParams *trackPar)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **fftFa1, **fftFa2;
	extern fftw_complex **psNoPad;
	extern fftw_complex **cNoPad;
	extern fftw_complex **cFast, **cFastOver;
	extern fftw_complex **img1, **img2;		 /* Detected images for amplitude match*/
	extern fftw_complex **img1in, **img2in;	/* ^^^ Detected images for amplitude match*/
	extern fftw_complex **fftFa1os, **fftFa2os; 
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **caNoPad;
	extern fftw_complex **psAmpNoPadL;
	extern float **caNoPadMag;
	extern double **meanS;
	extern double **sigmaS;
	extern double **corrResult;
	extern float **dataS, **dataR;
	fftw_complex **fftwTmp;
	strackComplex **ers1Tmp;
	int32_t wA2, wR2;
	int32_t i, j;
	/*
	  Input buffer.
	*/
	imageBuf1.nr = trackPar->imageP1.nSlpR;
	imageBuf1.na = NBUFFERLINES;
	imageBuf1.firstRow = LARGEINT;
	imageBuf1.lastRow = -1;
	imageBuf2.nr = trackPar->imageP2.nSlpR;
	imageBuf2.na = NBUFFERLINES;
	imageBuf2.firstRow = LARGEINT;
	imageBuf2.lastRow = -1;
	
	if (trackPar->floatFlag == TRUE)
	{
		fprintf(stderr, "Floating point input \n");
		imageBuf1.buf = (void **)mallocfftw_complexMat(imageBuf1.na, imageBuf1.nr);
		fftwTmp = (fftw_complex **)imageBuf1.buf;
		for (i = 0; i < NBUFFERLINES; i++)
			for (j = 0; j < imageBuf1.nr; j++)
			{
				fftwTmp[i][j].re = 1;
				fftwTmp[i][j].im = 0;
			}

		imageBuf2.buf = (void **)mallocfftw_complexMat(imageBuf2.na, imageBuf2.nr);
		fftwTmp = (fftw_complex **)imageBuf2.buf;
		for (i = 0; i < NBUFFERLINES; i++)
			for (j = 0; j < imageBuf2.nr; j++)
			{
				fftwTmp[i][j].re = 1;
				fftwTmp[i][j].im = 0;
			}
	}
	else
	{
		fprintf(stderr, "Short int input \n");
		imageBuf1.buf = (void *)mallocErs1ComplexMat(imageBuf1.na, imageBuf1.nr);
		ers1Tmp = (strackComplex **)imageBuf1.buf;
		for (i = 0; i < NBUFFERLINES; i++)
			for (j = 0; j < imageBuf1.nr; j++)
			{
				ers1Tmp[i][j].r = 1;
				ers1Tmp[i][j].i = 0;
			}

		imageBuf2.buf = (void *)mallocErs1ComplexMat(imageBuf2.na, imageBuf2.nr);
		ers1Tmp = (strackComplex **)imageBuf2.buf;
		for (i = 0; i < NBUFFERLINES; i++)
			for (j = 0; j < imageBuf2.nr; j++)
			{
				ers1Tmp[i][j].r = 1;
				ers1Tmp[i][j].i = 0;
			}
	}
	wA2 = (trackPar->wAa - 2 * trackPar->edgePadA) * OS;
	wR2 = (trackPar->wRa - 2 * trackPar->edgePadR) * OS;
	/*
	  Patches
	*/
	img1 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS); /* ^^^ */
	img2 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS); /* ^^^ */
	img1in = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);		  /* ^^^ */
	img2in = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);		  /* ^^^ */
	for (i = 0; i < trackPar->wAa * OS; i++)
		for (j = 0; j < trackPar->wRa * OS; j++)
		{
			img1[i][j].im = 0.;
			img2[i][j].im = 0.;
		}
	/*
	  FFTs
	*/
	fftFa1 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			   /* ^^^ */
	fftFa2 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			   /* ^^^ */
	fftFa1os = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);
	fftFa2os = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);
	psNoPad = mallocfftw_complexMat(trackPar->wA * OS, trackPar->wR * OS);
	psFast = mallocfftw_complexMat(NFAST+1, NFAST+1);
	psFastOver = mallocfftw_complexMat((NFAST+1) * NOVER, (NFAST+1) * NOVER);
	psAmpNoPad = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);
	for (i = 0; i < (NFAST+1) * NOVER; i++)
		for (j = 0; j < (NFAST+1) * NOVER; j++)
		{
			psFastOver[i][j].re = 0.0;
			psFastOver[i][j].im = 0.0;
		}
	cFast = mallocfftw_complexMat(NFAST + 1, NFAST + 1);
	cFastOver = mallocfftw_complexMat((NFAST + 1) * NOVER, (NFAST+1) * NOVER);
	caNoPad = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);
	caNoPadMag = mallocFloatMat(trackPar->wAa * OS, trackPar->wRa * OS);

	for (i = 0; i < (NFAST+1) * NOVER; i++)
		for (j = 0; j < (NFAST+1) * NOVER; j++)
		{
			cFastOver[i][j].re = 0.0;
			cFastOver[i][j].im = 0.0;
		}
	/* Input buffers */
	trackPar->offR = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->offA = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->corr = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->type = mallocByteMat(trackPar->nA, trackPar->nR);
	// clear buffers before starting 
	initMatrix(trackPar->offR, trackPar->nA, trackPar->nR, (float)(-LARGEINT));
	initMatrix(trackPar->offA, trackPar->nA, trackPar->nR, (float)(-LARGEINT));
	initMatrix(trackPar->corr, trackPar->nA, trackPar->nR, 0.0);
	initByteMatrix(trackPar->type, trackPar->nA, trackPar->nR, BAD);
	meanS = mallocDoubleMat(2 * trackPar->edgePadA * OS + 1, 2 * trackPar->edgePadR * OS + 1);
	sigmaS = mallocDoubleMat(2 * trackPar->edgePadA * OS + 1, 2 * trackPar->edgePadR * OS + 1);
	corrResult = mallocDoubleMat(2 * trackPar->edgePadA * OS + 1, 2 * trackPar->edgePadR * OS + 1);
	dataS = mallocFloatMat(trackPar->wAa * OS, trackPar->wRa * OS);
	initMatrix(dataS, trackPar->wAa * OS, trackPar->wRa * OS, 0.);
	dataR = mallocFloatMat(wA2, wR2);
	initMatrix(dataR, wA2, wR2, 0.);
}

/*************************************************************
   Output .dat file for matches
*************************************************************/
static void writeDatFile(TrackParams *trackPar)
{
	FILE *fp;
	fp = fopen(trackPar->outFileD, "w");
	fprintf(fp, "%i %i %i %i %i %i\n", trackPar->rStart * trackPar->scaleFactor, trackPar->aStart * trackPar->scaleFactor,
			trackPar->nR, trackPar->nA, trackPar->deltaR * trackPar->scaleFactor, trackPar->deltaA * trackPar->scaleFactor);
}

/*************************************************************
   Read mask
*************************************************************/
static int32_t maskValue(TrackParams *trackPar, int32_t r1, int32_t a1)
{
	int32_t ia1, ir1;
	int32_t dr, da;
	/* not mask, always try match */
	if (trackPar->maskFlag == FALSE)
		return (1);
	/* otherwise get mask value */
	/* fixed coords 10/29/2018 */
	dr = (trackPar->wRa / 2 - trackPar->maskDat.r0) * trackPar->scaleFactor;
	da = (trackPar->wAa / 2 - trackPar->maskDat.a0) * trackPar->scaleFactor;
	ia1 = (a1 + da) / trackPar->maskDat.nal;
	ir1 = (r1 + dr) / trackPar->maskDat.nrl;
	/* 	ia1 = (a1+trackPar->wA/2)/trackPar->maskDat.nal;	ir1 = (r1+trackPar->wR/2)/trackPar->maskDat.nrl; */
	if (ia1 < 0 || ia1 >= trackPar->maskDat.na || ir1 < 0 || ir1 >= trackPar->maskDat.nr)
	{
		return (0);
	}
	return ((int)(trackPar->maskDat.mask[ia1][ir1]));
}

/*************************************************************
   Routines to Malloc matrices
*************************************************************/
static float **mallocFloatMat(int32_t nA, int32_t nR)
{
	float *tmp, **tmp1;
	int32_t i;
	tmp = malloc(nR * nA * sizeof(float));
	tmp1 = (float **)malloc(nA * sizeof(float *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}

static double **mallocDoubleMat(int32_t nA, int32_t nR)
{
	double *tmp, **tmp1;
	int32_t i;
	tmp = malloc(nR * nA * sizeof(double));
	tmp1 = (double **)malloc(nA * sizeof(double *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
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

static fftw_complex **mallocfftw_complexMat(int32_t nA, int32_t nR)
{
	fftw_complex *tmp, **tmp1;
	int32_t i;
	tmp = malloc(nA * nR * sizeof(fftw_complex));
	tmp1 = (fftw_complex **)malloc(nA * sizeof(fftw_complex *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}

static strackComplex **mallocErs1ComplexMat(int32_t nA, int32_t nR)
{
	strackComplex *tmp, **tmp1;
	int32_t i;
	tmp = malloc(nA * nR * sizeof(strackComplex));
	tmp1 = (strackComplex **)malloc(nA * sizeof(strackComplex *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}

static void initMatrix(float **data, int32_t wA, int32_t wR, float value)
{
	int32_t i, j;
	for (i = 0; i < wA; i++)
		for (j = 0; j < wR; j++)
			data[i][j] = value;
}

static void initByteMatrix(char **data,  int32_t wA, int32_t wR, char value)
{
	int32_t i, j;
	for (i = 0; i < wA; i++)
		for (j = 0; j < wR; j++)
			data[i][j] = value;
}
