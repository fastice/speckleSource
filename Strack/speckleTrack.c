#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
#include <sys/types.h>
#include "time.h"
#include <libgen.h>
#include "cRecipes/nrutil.h"

/*#define NOVER 10
#define NFAST 16
*/
#define NOVER 10
#define NFAST 8
#define INTPATCHSIZE 8
#define LA 2
#define BAD 0
#define CMATCH 1
#define AMPMATCH 2
#define AMPMATCHLARGE 3
#define NBUFFERLINES 1000
#define FFTWPLANMODE FFTW_ESTIMATE
#define OS 2 /* Amplitude oversample factor */
#define OSA 2

static void ampMatch(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, int32_t large);
static void ampTrackFast(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, double p1, double p2, int32_t large);
static void cmpPSWithPad(TrackParams *trackPar);
static void cmpTrackFast(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax);
static void finalizeCmpxMatch(double cMax, int32_t iMax, int32_t jMax, double wRover2exp, double wAover2exp,
							  double invNOVER, TrackParams *trackPar, int32_t *good, int32_t *type,
							  float *rShift1, float *aShift1, double *cAvg, int32_t *nCorr);
static void findCPeak(float **corrF, fftw_complex **corrNoPad, int32_t wR, int32_t wA, int32_t *iMax, int32_t *jMax, double *cMax);
static void findImage2Pos(int32_t r1, int32_t a1, TrackParams *trackPar, int32_t *r2, int32_t *a2);
static void overSampleC(TrackParams *trackPar);
static void phaseCorrectBaseline(TrackParams *trackPar, int32_t r1, int32_t a1);
static void phaseCorrectInt(TrackParams *trackPar, int32_t r1, int32_t a1);
/* IO Routines */
static int32_t getAmpPatches(int32_t r1, int32_t a1, int32_t r2, int32_t a2, FILE *fp1, FILE *fp2, TrackParams *trackPar, int32_t large);
static int32_t getPatches(int32_t r1, int32_t a1, int32_t r2, int32_t a2, FILE *fp1, FILE *fp2, TrackParams *trackPar);
static void loadBuffer(StrackBuf *imageBuf1, SARData *imageP1, TrackParams *trackPar, int32_t a1, int32_t wAa, size_t sSize, FILE *fp1);
static void openStrackFiles(FILE **fpR, FILE **fpA, FILE **fpC, FILE **fpT, FILE **fp1, FILE **fp2, TrackParams *trackPar);
void printLineSummary(int32_t nTot, int32_t a1, double cAvg, double cAvgAmp, int32_t nMask, time_t lastTime,
					  time_t startTime, TrackParams *trackPar);
static void writeDatFile(TrackParams *trackPar);
/* Setup code */
static void clearTrackParBuffs(TrackParams *trackPar);
static void computeSinShift(double **sinShift, double **cosShift, TrackParams *trackPar);
static void fftPlans(TrackParams *trackPar);
static void initPhaseCorrect(TrackParams *trackPar);
static void mallocSpace(TrackParams *trackPar);
static void trackParInits(TrackParams *trackPar);
/* Utilities */
static float absCmpx(fftw_complex x, int32_t noSquareRoot);
static int32_t boundsCheck(int32_t r1, int32_t a1, int32_t r2, int32_t a2, int32_t wR, int32_t wA, TrackParams *trackPar);
static void computeHanning(TrackParams *trackPar);
static void getShifts(double range, double azimuth, Offsets *offsets, double *dr, double *da);
static int32_t inBounds(int32_t r1, int32_t a1, int32_t wR, int32_t wA, int32_t sR, int32_t sA);
static int32_t maskValue(TrackParams *trackPar, int32_t r1, int32_t a1);
static char **mallocByteMat(int32_t nA, int32_t nR);
static float **mallocFloatMat(int32_t nx, int32_t ny);
static fftw_complex **mallocfftw_complexMat(int32_t nx, int32_t ny);
static strackComplex **mallocErs1ComplexMat(int32_t nA, int32_t nR);
static void zeroComplex(fftw_complex **fft, int32_t na, int32_t nr);
static void zeroFloat(float **x, int32_t na, int32_t nr);
static void zeroPadFFT(fftw_complex **fftPad, int32_t pSA, int32_t pSR, fftw_complex **fftOrig, int32_t oSA, int32_t oSR, int32_t dA, int32_t dR);
/* Legacy Code */
static void cmpPSNoPad(TrackParams *trackPar);
static void estCarrier(TrackParams *trackPar);
static void estRangeCarrier(TrackParams *trackPar);
static void azComputeAzShift(int32_t j, int32_t *prevShift, TrackParams *trackPar, double *cosShift, double *sinShift);
/* Experimental code */
static void cmpTrackGauss(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax);
static void myGauss(xyData x, float a[], float *y, float dyda[], int32_t na);
void mrqcofMod(xyData x[], float y[], float sig[], int32_t ndata, float a[], int32_t ia[],
			   int32_t ma, float **alpha, float beta[], float *chisq,
			   void (*funcs)(xyData, float[], float *, float[], int));
void mrqminMod(xyData x[], float y[], float sig[], int32_t ndata, float a[], int32_t ia[],
			   int32_t ma, float **covar, float **alpha, float *chisq,
			   void (*funcs)(xyData, float[], float *, float[], int), float *alamda);

/* ************************ GLOBAL VARIABLES ***************************/
fftw_complex **fftF1, **fftF2; /* FFt's */
fftw_complex **patch1, **patch2;
fftw_complex **patch1in, **patch2in; /* ... */
fftw_complex **psNoPad;
fftw_complex **cFast, **cFastOver;
fftw_complex **cNoPad;
fftw_complex **psFast, **psFastOver;
fftw_complex **f1, **f2;
float **c, **cFastOverMag, **cNoPadMag;
float **hanning; /* Hanning window */
/* Amplitude stuff */
fftw_complex **img1, **img2; /* Detected images for amplitude match*/
fftw_complex **psAmpNoPad;
fftw_complex **fftFa1, **fftFa2; /* FFt's */
fftw_complex **caNoPad;
float **caNoPadMag;
/* Double Amplitude stuff */
fftw_complex **img1L, **img2L; /* Detected images for amplitude match*/
fftw_complex **psAmpNoPadL;
fftw_complex **fftFa1L, **fftFa2L; /* FFt's */
fftw_complex **caNoPadL;
float **caNoPadMagL;
/* Phase corrections */
double *thetaD;
double *sinThetaD;
double *cosThetaD;
double *overRange;
fftw_complex **intPatch;
fftw_complex **fftIntPatch;
fftw_complex **intPatchOver;
fftw_complex **fftIntPatchOver;
/* FFT plans */
fftwnd_plan aForward;
fftwnd_plan aForwardL;
fftwnd_plan aReverseNoPad;
fftwnd_plan aReverseNoPadL;
StrackBuf imageBuf1;
StrackBuf imageBuf2;
/* New stuff for oversampling ^^^ */
fftw_complex **img1in, **img2in;   /* ^^^ Detected images for amplitude match*/
fftw_complex **img1Lin, **img2Lin; /* ^^^ images for amplitude match*/
fftwnd_plan aForwardIn;
fftwnd_plan aForwardInL;
fftw_complex **fftFa1os, **fftFa1Los, **fftFa2os, **fftFa2Los; /*  ^^^ */
fftwnd_plan cForwardIn;										   /* ... */
fftw_complex **fftF1os, **fftF2os;							   /* ... */
/* ********************************************************************/
void speckleTrack(TrackParams *trackPar)
{
	time_t lastTime, myStart;
	FILE *fp1, *fp2, *fpR, *fpA, *fpC, *fpT;
	double *cosShift, *sinShift;
	double cMax, cAvg, cAvgAmp;
	double wRover2exp, wAover2exp, invNOVER, cThresh;
	float rShift, aShift, rShift1, aShift1;
	int32_t nTot, nCorr, nSkip, patchFlag, ampType;
	int32_t good, type, haveData, nMask, *nAmp, nAmpLine;
	int32_t i, j, r1, a1, r2, a2;
	int32_t r1a, a1a, r2a, a2a, jMax, iMax;
	int32_t prevShift, nTries, ampWScale, maskVal, large; /* Flag to determine whether to match */
	fprintf(stderr, "\nSPECKLE TRACKING\n");
	if (trackPar->legacyFlag == TRUE)
		trackPar->osF = OSA;
	else
		trackPar->osF = 1; /* Adjust cmpx over sample depending on mode */
	fprintf(stderr, "%li %li \n", (long) sizeof(fftw_real), (long) sizeof(fftw_complex));
	/*     Malloc space 	*/
	mallocSpace(trackPar);
	/* OPen files */
	openStrackFiles(&fpR, &fpA, &fpC, &fpT, &fp1, &fp2, trackPar);
	/* Init sin/cos shift ??? */
	computeSinShift(&sinShift, &cosShift, trackPar);
	/*  Compute plan for ffts and hanning window */
	fftPlans(trackPar);
	if (trackPar->hanningFlag == TRUE && trackPar->noComplex == FALSE)
		computeHanning(trackPar);
	/*  Initializations	*/
	wRover2exp = trackPar->wR / 2. * NOVER * OSA;
	wAover2exp = trackPar->wA / 2. * NOVER * OSA;
	invNOVER = 1.0 / (float)NOVER;
	trackParInits(trackPar);
	clearTrackParBuffs(trackPar); /* clear buffers before starting */
	/*
	  Loop to do matching
	*/
	prevShift = -1;
	nMask = 0;
	nTot = 0;
	trackPar->azShift = 0;
	myStart = time(NULL);
	for (i = 0; i < trackPar->nA; i++)
	{
		lastTime = time(NULL);
		a1 = trackPar->aStart + i * trackPar->deltaA - trackPar->wA / 2;
		nCorr = 0;
		nSkip = 0;
		cAvg = 0.0;
		cAvgAmp = 0.0;
		nAmpLine = 0;
		for (j = 0; j < trackPar->nR; j++)
		{
			haveData = TRUE;
			cMax = 0.;
			type = BAD;
			r1 = trackPar->rStart + j * trackPar->deltaR - trackPar->wR / 2;
			/* Find position for second image */
			findImage2Pos(r1, a1, trackPar, &r2, &a2);
			rShift = (r2 - r1);
			aShift = (a2 - a1);
			
			/* returns 1 if no mask available, 0 or 1 with mask */
			maskVal = maskValue(trackPar, r1 * trackPar->scaleFactor, a1 * trackPar->scaleFactor);
			if (maskVal > 1)
				maskVal = 1; /* This is because I may put larger values in mask for runcull */
			nMask += maskVal;
			/*
			   Complex matching
			*/
			if (boundsCheck(r1, a1, r2, a2, trackPar->wR, trackPar->wA, trackPar) && maskVal > 0 && trackPar->noComplex == FALSE)
			{
				
				/* Read patches from image */
				haveData = getPatches(r1, a1, r2, a2, fp1, fp2, trackPar);
				/*  Estimate any doppler component in azimuth */
				if (trackPar->legacyFlag == TRUE)
					azComputeAzShift(j, &prevShift, trackPar, cosShift, sinShift);
				good = FALSE;
				if (trackPar->noComplex == FALSE && haveData == TRUE)
				{
					phaseCorrectBaseline(trackPar, r1, a1); /* Phase correction with baseline  if present */
					phaseCorrectInt(trackPar, r1, a1);		/* Phase correction with interferogram if present */
					/* Compute initial cross correlation */
					if (trackPar->legacyFlag == TRUE)
						cmpPSNoPad(trackPar);
					else
						cmpPSWithPad(trackPar);
					/* Oversample peak; with oversamping or gaussian */
					if (trackPar->gaussFlag == TRUE)
						cmpTrackGauss(trackPar, &iMax, &jMax, &cMax);
					else
						cmpTrackFast(trackPar, &iMax, &jMax, &cMax);
					good = -1;
					finalizeCmpxMatch(cMax, iMax, jMax, wRover2exp, wAover2exp, invNOVER,
									  trackPar, &good, &type, &rShift1, &aShift1, &cAvg, &nCorr);
					
				}
				else
				{
					if (haveData == FALSE)
						nSkip++;
				}
			}
			else
				good = FALSE;
			/*
			  Amplitude matching
			*/
			nTries = 0; /* Counter for diff amp windows, should not exceed 2 */
			while (good == FALSE && haveData == TRUE && maskVal > 0 && nTries < trackPar->maxTries)
			{
				if (nTries == 0)
				{
					large = FALSE;
					ampWScale = 1;
					ampType = AMPMATCH;
					nAmp = &(trackPar->nAmp);
					cThresh = 0.07;
				}
				else
				{
					large = TRUE;
					ampWScale = LA;
					ampType = AMPMATCHLARGE;
					nAmp = &(trackPar->nAmpL);
					cThresh = 0.028;
				}
				a1a = trackPar->aStart + i * trackPar->deltaA - trackPar->wAa / 2 * ampWScale;
				r1a = trackPar->rStart + j * trackPar->deltaR - trackPar->wRa / 2 * ampWScale;
				r2a = r1a + (r2 - r1);
				a2a = a1a + (a2 - a1);
				
				if (boundsCheck(r1a, a1a, r2a, a2a, trackPar->wRa * ampWScale, trackPar->wAa * ampWScale, trackPar) && maskVal > 0)
				{
					haveData = getAmpPatches(r1a, a1a, r2a, a2a, fp1, fp2, trackPar, large);
					if (haveData == TRUE)
					{
						ampMatch(trackPar, &iMax, &jMax, &cMax, large);
						rShift1 = ((float)jMax - trackPar->wRa * 0.5 * OS * NOVER * ampWScale) * invNOVER / OS;
						aShift1 = ((float)iMax - trackPar->wAa * 0.5 * OS * NOVER * ampWScale) * invNOVER / OS;
					}
					else
					{
						cMax = 0.0;
						nSkip++;
					}
				}
				else
					cMax = 0.0;
				if (fabs((double)rShift1) < 0.15 * trackPar->wRa * OS * ampWScale &&
					fabs((double)aShift1) < 0.15 * trackPar->wAa * OS * ampWScale && cMax > cThresh)
				{
					good = TRUE;
					(*nAmp)++;
					type = ampType;
					cAvgAmp += cMax;
					nAmpLine++;
				}
				nTries++;
			}
			if (good == TRUE && maskVal > 0)
			{
				rShift -= rShift1;
				aShift -= aShift1;
			}
			else
			{
				rShift = -LARGEINT;
				aShift = -LARGEINT;
				if (haveData == TRUE)
					trackPar->nFail++;
			}
			/* Save values */
			trackPar->offR[i][j] = rShift;
			trackPar->offA[i][j] = aShift;
			if (rShift > -100000)
			{
				trackPar->offR[i][j] *= (float)trackPar->scaleFactor;
				trackPar->offA[i][j] *= (float)trackPar->scaleFactor;
			}
			trackPar->corr[i][j] = cMax; /*cMax;*/
			trackPar->type[i][j] = type; /*cMax;*/
		}								 /* End for j=0... */
		writeOffsets(i, trackPar, fpR, fpA, fpC, fpT);
		nTot += (trackPar->nR - nSkip);
		if (nCorr > 0)
		{
			cAvg = cAvg / (double)(nCorr);
		}
		else
			cAvg = 0.0;
		if (nAmpLine > 0)
			cAvgAmp /= nAmpLine;
		else
			cAvgAmp = 0.0;
		printLineSummary(nTot, a1, cAvg, cAvgAmp, nMask, lastTime, myStart, trackPar);
		//error("STOP");
	} /* End for i=0... */
	fprintf(stderr, "\n\nTotal time %f\n", (time(NULL) - myStart) / 60.);
	fclose(fpR);
	fclose(fpA);
	fclose(fpC);
	fclose(fpT);
	writeDatFile(trackPar);
	writeVrtFile(trackPar);
}

/**************************** Matching Code ***************************************/

/*
   Main part of amplitude matching.
 */
static void ampMatch(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, int32_t large)
{
	extern fftwnd_plan aForward;
	extern fftwnd_plan aForwardL;
	extern fftw_complex **fftFa1, **fftFa2;
	extern fftw_complex **fftFa1L, **fftFa2L;
	extern fftw_complex **img1, **img2;
	extern fftw_complex **img1L, **img2L;
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **psAmpNoPadL;
	fftwnd_plan *aFor;
	fftw_complex **im1, **im2, **fa1, **fa2;
	fftw_complex **psANoPad;
	double avgAmp1, avgAmp2, p1, p2;
	int32_t i, j, i1, wR, wA;
	if (large == TRUE)
	{
		fa1 = fftFa1L;
		fa2 = fftFa2L;
		im1 = img1L;
		im2 = img2L;
		psANoPad = psAmpNoPadL;
		wA = trackPar->wAa * LA;
		wR = trackPar->wRa * LA;
		aFor = &aForwardL;
		psANoPad = psAmpNoPadL;
	}
	else
	{
		fa1 = fftFa1;
		fa2 = fftFa2;
		im1 = img1;
		im2 = img2;
		psANoPad = psAmpNoPad;
		wA = trackPar->wAa;
		wR = trackPar->wRa;
		aFor = &aForward;
		psANoPad = psAmpNoPad;
	}
	/*   Step 1: Compute amp and subtract mean	*/
	avgAmp1 = 0.0;
	avgAmp2 = 0.0;
	for (i = 0; i < wA * OS; i++)
	{
		for (j = 0; j < wR * OS; j++)
		{
			im1[i][j].re = sqrt(im1[i][j].re);
			im2[i][j].re = sqrt(im2[i][j].re);
			avgAmp1 += im1[i][j].re;
			avgAmp2 += im2[i][j].re;
		}
	}
	/*  Step 2: Subtract mean and est avg pwr	*/
	avgAmp1 = avgAmp1 / (double)(wA * OS * wR * OS);
	avgAmp2 = avgAmp2 / (double)(wA * OS * wR * OS);
	p1 = 0.0;
	p2 = 0.0;
	for (i = 0; i < wA * OS; i++)
	{
		for (j = 0; j < wR * OS; j++)
		{
			if (im1[i][j].re > 1.0e-20)
				im1[i][j].re -= avgAmp1;
			im2[i][j].re -= avgAmp2;
			p1 += im1[i][j].re * im1[i][j].re;
			p2 += im2[i][j].re * im2[i][j].re;
		}
	}
	p1 = p1 / (double)(wA * OS * wR * OS);
	p2 = p2 / (wA * OS * wR * OS);
	/*  Step 3: FFT and power spectrum */
	fftwnd_one(*aFor, im1[0], fa1[0]);
	fftwnd_one(*aFor, im2[0], fa2[0]);
	/*  Step 4:   Zero pad fft for over sampling */
	i1 = 0;
	for (i = 0; i < wA * OS; i++)
	{
		for (j = 0; j < wR * OS; j++)
		{
			psANoPad[i][j].re = fa1[i1][j].re * fa2[i1][j].re + fa1[i1][j].im * fa2[i1][j].im;
			psANoPad[i][j].im = fa1[i1][j].im * fa2[i1][j].re - fa1[i1][j].re * fa2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
	/*  Step 5 Compute correlation fuction and find peak	*/
	ampTrackFast(trackPar, iMax, jMax, cMax, p1, p2, large);
}

/*
   Compute Correlation
 */
static void ampTrackFast(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax, double p1, double p2, int32_t large)
{
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **psAmpNoPadL;
	extern fftwnd_plan aReverseNoPad;
	extern fftwnd_plan aReverseNoPadL;
	extern fftw_complex **cFast, **cFastOver;
	extern fftw_complex **caNoPad;
	extern fftw_complex **caNoPadL;
	extern float **caNoPadMag;
	extern float **caNoPadMagL;
	fftwnd_plan *aReverseNP;
	fftw_complex **psANoPad;
	fftw_complex **caNoP;
	FILE *fp1;
	float **caNoPadM;
	double corr;
	int32_t i, j, i1, j1, iMax1, jMax1;
	int32_t wA, wR;
	if (large == TRUE)
	{
		psANoPad = psAmpNoPadL;
		wA = trackPar->wAa * LA;
		wR = trackPar->wRa * LA;
		aReverseNP = &aReverseNoPadL;
		psANoPad = psAmpNoPadL;
		caNoPadM = caNoPadMagL;
		caNoP = caNoPadL;
	}
	else
	{
		psANoPad = psAmpNoPad;
		wA = trackPar->wAa;
		wR = trackPar->wRa;
		aReverseNP = &aReverseNoPad;
		psANoPad = psAmpNoPad;
		caNoPadM = caNoPadMag;
		caNoP = caNoPad;
	}
	/* FFT to get correlation function	*/
	fftwnd_one(*aReverseNP, psANoPad[0], caNoP[0]);
	/* Find the peak */
	findCPeak(caNoPadM, caNoP, wR * OS, wA * OS, iMax, jMax, cMax);
	/* Return if no peak */
	if (*iMax < -999)
		return;
	/*  load from around peak to oversample;  */
	for (i = 0; i < NFAST; i++)
	{
		i1 = *iMax - NFAST / 2 + i;
		j1 = *jMax - NFAST / 2;
		for (j = 0; j < NFAST; j++)
		{
			cFast[i][j].re = caNoPadM[i1][j1];
			cFast[i][j].im = 0.0;
			j1++;
		}
	}
	/* Oversample data */
	overSampleC(trackPar);
	/*   Find max */
	iMax1 = 0;
	jMax1 = 0;
	*cMax = 0.0;
	for (i = (NFAST * NOVER) / 2 - NOVER * 2; i < (NFAST * NOVER) / 2 + NOVER * 2; i++)
		for (j = (NFAST * NOVER) / 2 - NOVER * 2; j < (NFAST * NOVER) / 2 + NOVER * 2; j++)
		{
			corr = (cFastOver[i][j].re * cFastOver[i][j].re + cFastOver[i][j].im * cFastOver[i][j].im);
			if (corr > *cMax)
			{
				*cMax = corr;
				iMax1 = i;
				jMax1 = j;
			}
		}
	/* Compute correlation */
	*cMax = sqrt(*cMax) / (pow((double)(wA * OS), 2) * pow((double)(wR * OS), 2) * pow((double)NFAST, 2));
	*cMax = *cMax / sqrt(p1 * p2);
	/* Compute offsets */
	*iMax = (*iMax * NOVER) + iMax1 - (NOVER * NFAST) / 2;
	*jMax = (*jMax * NOVER) + jMax1 - (NOVER * NFAST) / 2;
}

/*
Use non-oversampled complex data, then zero pad spectra to oversample correlation
*/
static void cmpPSWithPad(TrackParams *trackPar)
{
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **psNoPad;
	extern fftwnd_plan cForwardIn; /* ... */
	extern fftw_complex **fftF1os, **fftF2os;
	fftw_complex zero;
	int32_t i, j, i1, j1, i2, j2, ia, ja, j2a;
	/*	  Compute forward fft's	*/
	zero.re = 0.0;
	zero.im = 0.0;
	fftwnd_one(cForwardIn, patch1[0], fftF1os[0]);
	fftwnd_one(cForwardIn, patch2[0], fftF2os[0]);
	/* zero pad */
	for (i = 0; i < trackPar->wA * OSA; i++)
		for (j = 0; j < trackPar->wR * OSA; j++)
			psNoPad[i][j] = zero;
	for (i = 0; i < trackPar->wA / 2; i++)
	{
		i1 = OSA * trackPar->wA - trackPar->wA / OSA + i;
		i2 = (trackPar->wA / OSA + i + trackPar->azShift) % trackPar->wA;
		ia = (i + trackPar->azShift) % trackPar->wA;
		j1 = trackPar->wR * OSA - trackPar->wR / OSA;
		j2 = trackPar->wR / OSA + trackPar->rangeShift;
		for (j = 0; j < trackPar->wR / 2; j++)
		{
			ja = (j + trackPar->rangeShift) % trackPar->wR;
			j2a = j2 % trackPar->wR;
			psNoPad[i][j].re = fftF1os[ia][ja].re * fftF2os[ia][ja].re + fftF1os[ia][ja].im * fftF2os[ia][ja].im;
			psNoPad[i][j].im = fftF1os[ia][ja].im * fftF2os[ia][ja].re - fftF1os[ia][ja].re * fftF2os[ia][ja].im;
			psNoPad[i1][j].re = fftF1os[i2][ja].re * fftF2os[i2][ja].re + fftF1os[i2][ja].im * fftF2os[i2][ja].im;
			psNoPad[i1][j].im = fftF1os[i2][ja].im * fftF2os[i2][ja].re - fftF1os[i2][ja].re * fftF2os[i2][ja].im;
			psNoPad[i1][j1].re = fftF1os[i2][j2a].re * fftF2os[i2][j2a].re + fftF1os[i2][j2a].im * fftF2os[i2][j2a].im;
			psNoPad[i1][j1].im = fftF1os[i2][j2a].im * fftF2os[i2][j2a].re - fftF1os[i2][j2a].re * fftF2os[i2][j2a].im;
			psNoPad[i][j1].re = fftF1os[ia][j2a].re * fftF2os[ia][j2a].re + fftF1os[ia][j2a].im * fftF2os[ia][j2a].im;
			psNoPad[i][j1].im = fftF1os[ia][j2a].im * fftF2os[ia][j2a].re - fftF1os[ia][j2a].re * fftF2os[ia][j2a].im;
			j1++;
			j2++;
		} /* End for j */
	}	  /* End for i */
}

/*
  Directly oversample the area in a box of width NFAST around peak, with NOVER sample
 */
static void cmpTrackFast(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax)
{
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **cFast, **cFastOver;
	extern fftw_complex **cNoPad;
	extern float **cFastOverMag, **cNoPadMag;
	extern fftw_complex **patch1, **patch2;
	double hanningCorrection;
	double corr, osF;
	int32_t i, j, i1, j1, iMax1, jMax1;
	/* 1.5 removes bias introduced by single hanning window, added 10/19/18 */
	if (trackPar->hanningFlag == TRUE)
		hanningCorrection = 1.5;
	else
		hanningCorrection = 1.0;
	/*   FFT to get correlation function */
	fftwnd_one(trackPar->cReverseNoPad, psNoPad[0], cNoPad[0]);
	/*  Find correlation peak */
	*cMax = 0.;
	findCPeak(cNoPadMag, cNoPad, trackPar->wR * OSA, trackPar->wA * OSA, iMax, jMax, cMax);
	if (*iMax < -999)
		return;
	/*  Load from around peak to oversample;	*/
	for (i = 0; i < NFAST; i++)
	{
		i1 = *iMax - NFAST / 2 + i;
		j1 = *jMax - NFAST / 2;
		for (j = 0; j < NFAST; j++)
		{
			cFast[i][j].re = cNoPadMag[i1][j1];
			cFast[i][j].im = 0.0;
			j1++;
		}
	}
	/*  Oversample data*/
	overSampleC(trackPar);
	/*  Find max	*/
	iMax1 = 0;
	jMax1 = 0;
	*cMax = 0.0;
	for (i = (NFAST * NOVER) / 2 - NOVER * 2; i < (NFAST * NOVER) / 2 + NOVER * 2; i++)
		for (j = (NFAST * NOVER) / 2 - NOVER * 2; j < (NFAST * NOVER) / 2 + NOVER * 2; j++)
		{
			corr = (cFastOver[i][j].re * cFastOver[i][j].re + cFastOver[i][j].im * cFastOver[i][j].im);
			if (corr > *cMax)
			{
				*cMax = corr;
				iMax1 = i;
				jMax1 = j;
			}
		}
	/* Added 10/19/18 for hanning bias */
	if (trackPar->legacyFlag == TRUE)
		osF = OSA;
	else
		osF = 1.;
	*cMax = hanningCorrection * sqrt(*cMax) /
			(pow((double)(trackPar->wA * osF), 2) * pow((double)(trackPar->wR * osF), 2) * pow((double)NFAST, 2));
	*cMax = *cMax / sqrt(trackPar->p1 * trackPar->p2);
	*iMax = (*iMax * NOVER) + iMax1 - (NOVER * NFAST) / 2;
	*jMax = (*jMax * NOVER) + jMax1 - (NOVER * NFAST) / 2;
}

/*
  Do the final calculations to scale match and decide if its good
*/
static void finalizeCmpxMatch(double cMax, int32_t iMax, int32_t jMax, double wRover2exp, double wAover2exp,
							  double invNOVER, TrackParams *trackPar, int32_t *good, int32_t *type, float *rShift1, float *aShift1,
							  double *cAvg, int32_t *nCorr)
{
	/* Update stats */
	if (cMax > 0)
	{
		*cAvg += cMax;
		(*nCorr)++;
	}
	*good = FALSE;
	if (cMax > .185 && trackPar->noComplex == FALSE)
	{
		/* Update shift with correction */
		*rShift1 = ((float)jMax - wRover2exp) * invNOVER / OSA;
		*aShift1 = ((float)iMax - wAover2exp) * invNOVER / OSA;
		/* Check that mach is within 1/8 of the window size */
		if (fabs((double)*rShift1) < 0.125 * trackPar->wR * OSA &&
			fabs((double)*aShift1) < 0.125 * trackPar->wA * OSA)
		{
			*good = TRUE;
			trackPar->nComplex++;
			*type = CMATCH;
		}
	}
}

/*
  Find peak of no pad correlation function
*/
static void findCPeak(float **corrF, fftw_complex **corrNoPad,
					  int32_t wR, int32_t wA, int32_t *iMax, int32_t *jMax, double *cMax)
{
	int32_t i, j, i1, j1;
	/*	  Find max value	*/
	*cMax = 0.0;
	for (i = 0; i < wA / 2; i++)
	{
		i1 = (wA) / 2 + i;
		j1 = (wR) / 2;
		for (j = 0; j < wR / 2; j++)
		{
			corrF[i][j] = absCmpx(corrNoPad[i1][j1], FALSE);
			if (corrF[i][j] > *cMax)
			{
				*cMax = corrF[i][j];
				*iMax = i;
				*jMax = j;
			}

			corrF[i1][j1] = absCmpx(corrNoPad[i][j], FALSE);
			if (corrF[i1][j1] > *cMax)
			{
				*cMax = corrF[i1][j1];
				*iMax = i1;
				*jMax = j1;
			}

			corrF[i1][j] = absCmpx(corrNoPad[i][j1], FALSE);
			if (corrF[i1][j] > *cMax)
			{
				*cMax = corrF[i1][j];
				*iMax = i1;
				*jMax = j;
			}

			corrF[i][j1] = absCmpx(corrNoPad[i1][j], FALSE);
			if (corrF[i][j1] > *cMax)
			{
				*cMax = corrF[i][j1];
				*iMax = i;
				*jMax = j1;
			}
			j1++;
		}
	}
	/*  Check that peak is within range */
	if (*iMax < (NFAST / 2) || *iMax > (wA - NFAST / 2) || *jMax < (NFAST / 2) || *jMax > (wR - NFAST / 2))
	{
		*cMax = -1.0;
		*jMax = -1000.;
		*iMax = -1000.;
	}
}

/*
  Find inital shift either from polynomial or initial guess.
*/
static void findImage2Pos(int32_t r1, int32_t a1, TrackParams *trackPar, int32_t *r2, int32_t *a2)
{
	double aShift, rShift;
	double sc;
	int32_t haveShift;
	haveShift = FALSE;
	sc = (double)trackPar->scaleFactor;
	if (trackPar->polyShift == FALSE)
	{
		getShifts(r1, a1, &(trackPar->initOffsets), &rShift, &aShift);
		if (rShift > (-LARGEINT + 1))
		{
			haveShift = TRUE;
		}
	}
	if (trackPar->polyShift == TRUE || haveShift == FALSE)
	{
		rShift = trackPar->rShiftPoly[0] + (double)(r1 * sc) * trackPar->rShiftPoly[1] +
				 (double)(a1 * sc) * trackPar->rShiftPoly[2];
		aShift = trackPar->aShiftPoly[0] + (double)(r1 * sc) * trackPar->aShiftPoly[1] +
				 (double)(a1 * sc) * trackPar->aShiftPoly[2];
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

/*
   Over sample area around peak.
 */
static void overSampleC(TrackParams *trackPar)
{
	extern fftw_complex **psFast, **psFastOver;
	extern fftw_complex **cFast, **cFastOver;
	int32_t i, j;
	int32_t i1, i2, j1, j2;
	fftwnd_one(trackPar->cForwardFast, cFast[0], psFast[0]);
	/* Zero pad fft for over sampling	*/
	zeroPadFFT(psFastOver, NFAST * NOVER, NFAST * NOVER, psFast, NFAST, NFAST, 0, 0);
	/* Inverse transform */
	fftwnd_one(trackPar->cReverseFast, psFastOver[0], cFastOver[0]);
}

/*
Use baseline to estimate and remove phase ramp.
*/
static void phaseCorrectBaseline(TrackParams *trackPar, int32_t r1, int32_t a1)
{
	extern double *sinThetaD;
	extern double *cosThetaD;
	extern double *overRange;
	double x, deltaX, ptod, bnPhase;
	double pr, pi, cr, ci;
	double bn, bp, bsq;
	int32_t i, j, j1, j2, i1, i2;
	int32_t overSample;
	/* account for new patch oversampling */
	if (trackPar->legacyFlag == FALSE)
		overSample = 1;
	else
		overSample = OSA;
	/* Along track coordate (-.5,.5) */
	x = (a1 - trackPar->wA / 2 - (trackPar->imageP1.nSlpA / 2)) /
		(double)(trackPar->imageP1.nSlpA);
	deltaX = 1. / (double)(trackPar->imageP1.nSlpA * overSample);
	bn = trackPar->Bn + x * trackPar->dBn + x * x * trackPar->dBnQ;
	ptod = -(4.0 * PI) / trackPar->lambda;
	/*
	  Baseline correction
	*/
	for (i = 0; i < trackPar->wA * overSample; i++)
	{
		bp = trackPar->Bp + x * trackPar->dBp + x * x * trackPar->dBpQ;
		bsq = bn * bn + bp * bp;
		for (j = 0; j < trackPar->wR * overSample; j++)
		{
			j1 = max(r1 * overSample - trackPar->wR * overSample / 2, 0) + j;
			bnPhase = -bn * sinThetaD[j1] - bp * cosThetaD[j1] + (0.5 * bsq * overRange[j1]);
			bnPhase *= ptod;
			pr = cos(bnPhase);
			pi = sin(bnPhase);
			cr = patch1[i][j].re;
			ci = patch1[i][j].im;
			patch1[i][j].re = cr * pr - ci * pi;
			patch1[i][j].im = cr * pi + ci * pr;
		}
		x += deltaX;
	}
}

static void phaseCorrectInt(TrackParams *trackPar, int32_t r1, int32_t a1)
{
	extern fftw_complex **intPatch;
	extern fftw_complex **fftIntPatch;
	extern fftw_complex **intPatchOver;
	extern fftw_complex **fftIntPatchOver;
	extern fftw_complex **intPatch, **fftIntPatch;
	extern fftw_complex **patch1, **patch2;
	double cr, ci, pr, pi;
	double maxf;
	double tmp;
	double p, scale;
	int32_t i, j, j1, j2, i1, i2;
	int32_t patchSize;
	double lambda;
	int32_t ir1, ia1;
	patchSize = trackPar->intDat.patchSize;
	/*
	  Use interferogram to improve match
	*/
	if (trackPar->intFlag == TRUE)
	{
		ia1 = (a1 + trackPar->wA / 2) / trackPar->intDat.nal - trackPar->intDat.patchSize / 2;
		ir1 = (r1 + trackPar->wR / 2) / trackPar->intDat.nrl - trackPar->intDat.patchSize / 2;

		if (inBounds(ir1, ia1, patchSize, patchSize, trackPar->intDat.nr, trackPar->intDat.na))
		{
			/* Read patch of interferogram data*/
			for (i = 0; i < trackPar->intDat.patchSize; i++)
			{
				for (j = 0; j < trackPar->intDat.patchSize; j++)
				{
					intPatch[i][j].re = trackPar->intDat.intf[ia1 + i][ir1 + j].re;
					intPatch[i][j].im = trackPar->intDat.intf[ia1 + i][ir1 + j].im;
				}
			}
			/* FFT PATCH */
			fftwnd_one(trackPar->intDat.forward, intPatch[0], fftIntPatch[0]);
			/* Find Peak in interferogram */
			maxf = 0.;
			for (i = 0; i < trackPar->intDat.patchSize; i++)
			{
				for (j = 0; j < trackPar->intDat.patchSize; j++)
				{
					p = fftIntPatch[i][j].re * fftIntPatch[i][j].re + fftIntPatch[i][j].im * fftIntPatch[i][j].im;
					maxf = max(maxf, p);
				}
			}
			/* Zero out everything but the peak, which will have the fringe ramp */
			for (i = 0; i < trackPar->intDat.patchSize; i++)
			{
				for (j = 0; j < trackPar->intDat.patchSize; j++)
				{
					p = fftIntPatch[i][j].re * fftIntPatch[i][j].re + fftIntPatch[i][j].im * fftIntPatch[i][j].im;
					if (p < 0.99 * maxf)
						fftIntPatch[i][j].re = 0.0;
					fftIntPatch[i][j].im = 0.0;
				}
			}
			/* OVer Sampling via fft and zero pad */
			zeroPadFFT(fftIntPatchOver, patchSize * trackPar->intDat.nal * trackPar->osF,
					   patchSize * trackPar->intDat.nrl * trackPar->osF, fftIntPatch, patchSize, patchSize, 0, 0);
			/* Inverse transform */
			fftwnd_one(trackPar->intDat.backward, fftIntPatchOver[0], intPatchOver[0]);
			/* NOw do correction */
			i1 = (patchSize * trackPar->intDat.nal * trackPar->osF) / 2 - trackPar->wA * trackPar->osF / 2;
			for (i = 0; i < trackPar->wA * trackPar->osF; i++)
			{
				j1 = (patchSize * trackPar->intDat.nrl * trackPar->osF) / 2 - trackPar->wR * trackPar->osF / 2;
				for (j = 0; j < trackPar->wR * trackPar->osF; j++)
				{
					pr = intPatchOver[i1][j1].re;
					pi = -intPatchOver[i1][j1].im;
					if (fabs((double)pr) < 1.0e-6 && fabs((double)pi) < 1.0e-6)
					{
						pr = 1.0;
						pi = 0.0;
						scale = 1.0;
					}
					else
						scale = 1.0 / sqrt(pr * pr + pi * pi);
					cr = patch1[i][j].re;
					ci = patch1[i][j].im;
					patch1[i][j].re = (cr * pr - ci * pi) * scale;
					patch1[i][j].im = (cr * pi + ci * pr) * scale;
					j1++;
				}
				i1++;
			}
		} /* End if ia1 > 0 && */
	}	  /* End if(trackPar->intFlag == TRUE  */
}

/* ************************ IO Routines ********************************/
/*
  Read patches for amplitude matching
*/
static int32_t getAmpPatches(int32_t r1, int32_t a1, int32_t r2, int32_t a2, FILE *fp1, FILE *fp2, TrackParams *trackPar, int32_t large)
{
	size_t offset1, offset2;
	int32_t i, j, i1, j1, i2, j2, ja, j2a;
	int32_t dOff1, dOff2;
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	extern fftw_complex **img1, **img2;
	extern fftw_complex **img1L, **img2L;
	extern fftw_complex **img1in, **img2in;								  /* ^^^ Detected images for amplitude match*/
	extern fftw_complex **img1Lin, **img2Lin;							  /* ^^^ images for amplitude match*/
	extern fftwnd_plan aForwardIn;										  /* ^^^ */
	extern fftwnd_plan aForwardInL;										  /* ^^^ */
	extern fftwnd_plan aReverseNoPadL;									  /* ^^^ */
	extern fftwnd_plan aReverseNoPadL;									  /* ^^^ */
	extern fftw_complex **fftFa1os, **fftFa1Los, **fftFa2os, **fftFa2Los; /*  ^^^ */
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **fftFa1, **fftFa2;	  /* FFt's */
	extern fftw_complex **fftFa1L, **fftFa2L; /* FFt's */
	fftw_complex **im1, **im2, **im1in, **im2in;
	strackComplex **ers1Buf1, **ers1Buf2;
	fftw_complex **fftwBuf1, **fftwBuf2;
	fftw_complex **f1, **f2, **f1a, **f2a;
	fftw_complex zero;
	size_t sSize, wRa, wAa, s1, s2, a1a, a2a;
	double i1Sum, i2Sum;
	float scale;
	int32_t azShift, rangeShift;
	int32_t ia, i2a, nZero;
	void *tmp;
	zero.re = 0.0;
	zero.im = 0.0;
	nZero = 0;
	if (a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA)
	{
		error("getpatches");
		return (FALSE);
	}
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
	if (large == TRUE)
	{
		wAa = trackPar->wAa * LA;
		wRa = trackPar->wRa * LA;
		im1 = img1L;
		im2 = img2L;
		im1in = img1Lin;
		im2in = img2Lin;
		f1 = fftFa1Los;
		f2 = fftFa2Los;
		f1a = fftFa1L;
		f2a = fftFa2L;
	}
	else
	{
		wAa = trackPar->wAa;
		wRa = trackPar->wRa;
		im1 = img1;
		im2 = img2;
		im1in = img1in;
		im2in = img2in;
		f1 = fftFa1os;
		f2 = fftFa2os;
		f1a = fftFa1;
		f2a = fftFa2;
	}
	scale = 1. / (wRa * wAa * wRa * wAa * OS * OS * OS); /* rough guess at scale to avoid fp overflow */
	/*
	  Load buffers if needed
	*/
	loadBuffer(&imageBuf1, &(trackPar->imageP1), trackPar, a1, wAa, sSize, fp1);
	loadBuffer(&imageBuf2, &(trackPar->imageP2), trackPar, a2, wAa, sSize, fp2);
	/*
	  Get complex data
	*/
	for (i = 0; i < wAa; i++)
	{
		i1 = a1 - imageBuf1.firstRow + i;
		i2 = a2 - imageBuf2.firstRow + i;
		for (j = 0; j < wRa; j++)
		{
			j1 = r1 + j;
			j2 = r2 + j;
			if (trackPar->floatFlag == TRUE)
			{
				im1in[i][j].re = (fftw_real)fftwBuf1[i1][j1].re;
				im2in[i][j].re = (fftw_real)fftwBuf2[i2][j2].re;
				im1in[i][j].im = (fftw_real)fftwBuf1[i1][j1].im;
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
	/* Forward FFT */
	if (large == FALSE)
	{
		/*fprintf(stderr,"forward fft \n");*/
		fftwnd_one(aForwardIn, im1in[0], f1[0]);
		fftwnd_one(aForwardIn, im2in[0], f2[0]);
	}
	else
	{
		fftwnd_one(aForwardInL, im1in[0], f1[0]);
		fftwnd_one(aForwardInL, im2in[0], f2[0]);
	}
	azShift = (int)((float)trackPar->azShift * (float)wAa / (float)trackPar->wA);
	rangeShift = (int)((float)trackPar->rangeShift * (float)wRa / (float)trackPar->wR);
	/*
	  Zero pad fft for over sampling
	*/
	zeroComplex(f1a, wAa * OS, wRa * OS);
	zeroComplex(f2a, wAa * OS, wRa * OS);
	zeroPadFFT(f1a, wAa * OS, wRa * OS, f1, wAa, wRa, 0, 0);
	zeroPadFFT(f2a, wAa * OS, wRa * OS, f2, wAa, wRa, 0, 0);
	/* Inverse transform */
	if (large == FALSE)
	{
		fftwnd_one(aReverseNoPad, f1a[0], im1[0]);
		fftwnd_one(aReverseNoPad, f2a[0], im2[0]);
	}
	else
	{
		fftwnd_one(aReverseNoPadL, f1a[0], im1[0]);
		fftwnd_one(aReverseNoPadL, f2a[0], im2[0]);
	}
	/*	  Detect Data 	*/
	i1Sum = 0.0;
	i2Sum = 0.0;
	for (i = 0; i < wAa * OS; i++)
	{
		for (j = 0; j < wRa * OS; j++)
		{
			im1[i][j].re = absCmpx(im1[i][j], TRUE) * scale;
			im1[i][j].im = 0.0;
			im2[i][j].re = absCmpx(im2[i][j], TRUE) * scale;
			im2[i][j].im = 0.0;
			if (im1[i][j].re < 1.e-9 || im2[i][j].re < 1e-9)
				nZero++;
			i1Sum += im1[i][j].re;
			i2Sum += im2[i][j].re;
		}
	}
	if (nZero > 20)
		return (FALSE);
	/* Return false if either data is just zeros; added 8/11/2016 to avoid zeros in tops images */
	if (i1Sum < 1.0e-6 || i2Sum < 1.0e-6)
		return (FALSE);
	return (TRUE);
}

/*
  Read patches for complex matching
*/
static int32_t getPatches(int32_t r1, int32_t a1, int32_t r2, int32_t a2, FILE *fp1, FILE *fp2, TrackParams *trackPar)
{
	size_t offset1, offset2;
	size_t dOff1, dOff2;
	int32_t i, j, i1, j1, i2, j2, ia, ja, j2a;
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **patch1in, **patch2in;
	extern fftwnd_plan cForwardIn; /* ... */
	fftw_complex **patch1Use, **patch2Use;
	strackComplex **ers1Buf1, **ers1Buf2;
	fftw_complex **fftwBuf1, **fftwBuf2;
	fftw_complex zero;
	double p1, p2;
	float scale, power1, power2;
	int32_t nZero;
	size_t sSize, wR, s1, s2, a1a, a2a;
	zero.re = 0.0;
	zero.im = 0.0;
	nZero = 0;
	if (a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA)
	{
		error("getpatches");
		return (FALSE);
	}
	wR = trackPar->wR;
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
	/*
	  Load buffers if needed
	*/
	loadBuffer(&imageBuf1, &(trackPar->imageP1), trackPar, a1, trackPar->wA, sSize, fp1);
	loadBuffer(&imageBuf2, &(trackPar->imageP2), trackPar, a2, trackPar->wA, sSize, fp2);
	/*   Loop to input patches    */
	if (trackPar->legacyFlag == TRUE)
	{
		patch1Use = patch1in;
		patch2Use = patch2in;
	}
	else
	{
		patch1Use = patch1;
		patch2Use = patch2;
	}
	/*
	   Non-legacy case overSample in correlation not here.
	 */
	p1 = 0.;
	p2 = 0.0;
	for (i = 0; i < trackPar->wA; i++)
	{
		i1 = a1 - imageBuf1.firstRow + i;
		i2 = a2 - imageBuf2.firstRow + i;
		for (j = 0; j < trackPar->wR; j++)
		{
			j1 = r1 + j;
			j2 = r2 + j;
			if (trackPar->floatFlag == TRUE)
			{
				patch1Use[i][j].re = fftwBuf1[i1][j1].re * 1000;
				patch1Use[i][j].im = fftwBuf1[i1][j1].im * 1000;
				patch2Use[i][j].re = fftwBuf2[i2][j2].re * 1000;
				patch2Use[i][j].im = fftwBuf2[i2][j2].im * 1000;
			}
			else
			{
				patch1Use[i][j].re = (fftw_real)(ers1Buf1[i1][j1].r);
				patch2Use[i][j].re = (fftw_real)(ers1Buf2[i2][j2].r);
				patch1Use[i][j].im = (fftw_real)(ers1Buf1[i1][j1].i);
				patch2Use[i][j].im = (fftw_real)(ers1Buf2[i2][j2].i);
			}
			
			if (trackPar->hanningFlag == TRUE)
			{
				/* removed patch 1 weighting 10/19/18, compensated correlation with 1.5 scale */
				patch2Use[i][j].re *= hanning[i][j];
				patch2Use[i][j].im *= hanning[i][j];
			}
			power1 = absCmpx(patch1Use[i][j], TRUE);
			power2 = absCmpx(patch2Use[i][j], TRUE);
			if (power1 < 1e-9 || power2 < 1e-9) {
				nZero++;
			}
			
			/* If more than 20 effectively zero points, skip patch - avoids edges */
			if (nZero > 20) {
				return (FALSE);
			}
				
			p1 += power1;
			p2 += power2;
		}
	}
	if (p1 < 1e-6 || p2 < 1e-6)
		return (FALSE);
	if (trackPar->legacyFlag == FALSE)
	{
		/* Scale and compute power */
		trackPar->p1 = p1 / (double)(trackPar->wA * trackPar->wR);
		trackPar->p2 = p2 / (double)(trackPar->wA * trackPar->wR);
		return (TRUE);
	}
	/*
	   Legacy case. Oversample here and apply some other corrections.
	*/
	estCarrier(trackPar);
	estRangeCarrier(trackPar);
	/* Forward FFT */
	fftwnd_one(cForwardIn, patch1in[0], fftF1os[0]);
	fftwnd_one(cForwardIn, patch2in[0], fftF2os[0]);
	/* Zero arrays */
	zeroPadFFT(fftF1, OSA * trackPar->wA, OSA * trackPar->wR, fftF1os, trackPar->wA, trackPar->wR, trackPar->azShift, trackPar->rangeShift);
	zeroPadFFT(fftF2, OSA * trackPar->wA, OSA * trackPar->wR, fftF1os, trackPar->wA, trackPar->wR, trackPar->azShift, trackPar->rangeShift);
	/* Inverse fft */
	fftwnd_one(trackPar->cReverseNoPad, fftF1[0], patch1[0]);
	fftwnd_one(trackPar->cReverseNoPad, fftF2[0], patch2[0]);
	/* Scale and compute power */
	p1 = 0.;
	p2 = 0.0;
	scale = 1.0 / (trackPar->wA * trackPar->wR * OSA);
	for (i = 0; i < trackPar->wA * OSA; i++)
	{
		for (j = 0; j < trackPar->wR * OSA; j++)
		{
			patch1[i][j].re *= scale;
			patch1[i][j].im *= scale;
			patch2[i][j].re *= scale;
			patch2[i][j].im *= scale;
			p1 += absCmpx(patch1[i][j], TRUE);
			p2 += absCmpx(patch2[i][j], TRUE);
		}
	}
	trackPar->p1 = p1 / (double)(trackPar->wA * OSA * trackPar->wR * OSA);
	trackPar->p2 = p2 / (double)(trackPar->wA * OSA * trackPar->wR * OSA);
	return (TRUE);
}

/*
  Check if buffer needs refreshing.
 */
static void loadBuffer(StrackBuf *imageBuf1, SARData *imageP1, TrackParams *trackPar,
					   int32_t a1, int32_t wAa, size_t sSize, FILE *fp1)
{
	size_t a1a, s1, offset1;

	if (a1 < imageBuf1->firstRow || (a1 + wAa) > imageBuf1->lastRow)
	{
		a1a = max(0, a1 - trackPar->wAa * LA);
		s1 = imageP1->nSlpR * min(imageP1->nSlpA - a1a, NBUFFERLINES);
		offset1 = ((off_t)a1a * imageP1->nSlpR) * (off_t)sSize;
		fseeko(fp1, offset1, SEEK_SET);
		if (trackPar->floatFlag == TRUE) {
			if(trackPar->byteOrder == MSB)
			{
				freadBS((void *)imageBuf1->buf[0], sSize, s1, fp1, FLOAT32FLAG);
			} else {
				size_t rv = fread((void *)imageBuf1->buf[0], sSize, s1, fp1);
			}
		} else 
		{
			freadBS((void *)imageBuf1->buf[0], sSize, s1, fp1, INT16FLAG);	
		}		
		imageBuf1->firstRow = a1a;
		imageBuf1->lastRow = a1a + min(imageP1->nSlpA - a1a, NBUFFERLINES) - 1;
	}
}

/*
open all the necessary files
 */
static void openStrackFiles(FILE **fpR, FILE **fpA, FILE **fpC, FILE **fpT, FILE **fp1, FILE **fp2, TrackParams *trackPar)
{
	*fpR = fopen(trackPar->outFileR, "w");
	*fpA = fopen(trackPar->outFileA, "w");
	*fpC = fopen(trackPar->outFileC, "w");
	*fpT = fopen(trackPar->outFileT, "w");
	/*
		  Open input files
	*/
	*fp1 = openInputFile(trackPar->imageFile1);
	*fp2 = openInputFile(trackPar->imageFile2);
}

void printLineSummary(int32_t nTot, int32_t a1, double cAvg, double cAvgAmp, int32_t nMask, time_t lastTime,
					  time_t startTime, TrackParams *trackPar)
{
	int32_t totalMatches, totalAmp, lineTime;
	double percentMatch, percentComplex, percentAmp, percentFail, totalTime;

	totalMatches = trackPar->nComplex + trackPar->nAmp + trackPar->nAmpL;
	/* Percent matches of those attempted */
	percentMatch = 100. * (double)(trackPar->nComplex + trackPar->nAmp + trackPar->nAmpL) / (double)nTot;
	percentComplex = 100. * (double)(trackPar->nComplex) / (double)(nTot);
	totalAmp = trackPar->nAmp + trackPar->nAmpL;
	percentAmp = 100. * (double)(trackPar->nAmp + trackPar->nAmpL) / (double)(nTot);
	percentFail = 100. * (double)(trackPar->nFail) / (double)(nTot);
	lineTime = (int)(time(NULL) - lastTime);
	totalTime = (double)(time(NULL) - startTime) / 60.;
	fprintf(stderr,
			"\r%6i Tot %i, nMtch %8i %4.1f nCmpx %8i %4.1f "
			"nAmp %8i %4.1f (%7i -%7i) nFl %7i %4.1f cAvg %4.2f %4.2f %i t %3i(s)  %5.1f(m)",
			a1,
			nTot, totalMatches, percentMatch, trackPar->nComplex, percentComplex, totalAmp, percentAmp,
			trackPar->nAmp, trackPar->nAmpL, trackPar->nFail, percentFail, cAvg, cAvgAmp, nMask, lineTime, totalTime);
}


/*
   Write offsets data file
 */
static void writeDatFile(TrackParams *trackPar)
{
	FILE *fp;
	char geo2[2048], *tmp, filename[2048];
	fp = fopen(trackPar->outFileD, "w");
	fprintf(fp, "%i %i %i %i %i %i\n", trackPar->rStart * trackPar->scaleFactor, trackPar->aStart * trackPar->scaleFactor,
			trackPar->nR, trackPar->nA, trackPar->deltaR * trackPar->scaleFactor, trackPar->deltaA * trackPar->scaleFactor);
	// Make a copy of file name so dirname doesnot corrupt
	filename[0] = '\0';
	strcpy(filename, trackPar->imageFile2);
	/* Add geodats for pair */
	geo2[0] = '\0';
	tmp = strcat(geo2, dirname(filename));
	tmp = strcat(geo2, "/");
	tmp = strcat(geo2, trackPar->intGeodat);
	fprintf(fp, "%s %s\n", trackPar->intGeodat, geo2);
}


/* ************************ Setup Routines ********************************/
/*
   Initialize buffers
*/
static void clearTrackParBuffs(TrackParams *trackPar)
{
	int32_t i, j;
	for (i = 0; i < trackPar->nA; i++)
	{
		for (j = 0; j < trackPar->nR; j++)
		{
			trackPar->offR[i][j] = (float)(-LARGEINT);
			trackPar->offA[i][j] = (float)(-LARGEINT);
			trackPar->corr[i][j] = 0;
			trackPar->type[i][j] = BAD;
		}
	}
}


/*
  Pre compute sinShitfs. (legacy)
 */
static void computeSinShift(double **sinShift, double **cosShift, TrackParams *trackPar)
{
	int32_t i;
	double tmpx1;
	*sinShift = malloc(trackPar->wA * sizeof(double));
	*cosShift = malloc(trackPar->wA * sizeof(double));
	for (i = 0; i < trackPar->wA; i++)
	{
		tmpx1 = 2.0 * PI / (double)(trackPar->wA) * (double)i;
		(*sinShift)[i] = sin(tmpx1);
		(*cosShift)[i] = cos(tmpx1);
	}
}

/*
  Do plans for fft's
*/
static void fftPlans(TrackParams *trackPar)
{
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **f1, **f2;
	extern fftwnd_plan aForward;
	extern fftwnd_plan aForwardL;
	extern fftwnd_plan aForwardIn;
	extern fftwnd_plan cForwardIn;
	extern fftwnd_plan aForwardInL;
	extern fftwnd_plan aReverseNoPad;
	extern fftwnd_plan aReverseNoPadL;

	trackPar->cForward = fftw2d_create_plan(trackPar->wA * OSA, trackPar->wR * OSA,
											FFTW_FORWARD, FFTWPLANMODE);
	cForwardIn = fftw2d_create_plan(trackPar->wA, trackPar->wR,
									FFTW_FORWARD, FFTWPLANMODE);
	fprintf(stderr, "1 %i %i\n", NFAST, NFAST);

	trackPar->cForwardFast = fftw2d_create_plan(NFAST, NFAST,
												FFTW_FORWARD, FFTWPLANMODE);
	fprintf(stderr, "2 %i %i\n", trackPar->wA * OSA, trackPar->wR * OSA);
	trackPar->cReverseNoPad = fftw2d_create_plan(trackPar->wA * OSA,
												 trackPar->wR * OSA, FFTW_BACKWARD, FFTWPLANMODE);
	fprintf(stderr, "3 %i %i\n", NOVER * NFAST, NOVER * NFAST);
	trackPar->cReverseFast = fftw2d_create_plan(NOVER * NFAST,
												NOVER * NFAST, FFTW_BACKWARD, FFTWPLANMODE);
	fprintf(stderr, "4 %i\n", trackPar->wA);
	trackPar->onedForward1 = fftw_create_plan_specific(trackPar->wA,
													   FFTW_FORWARD, FFTWPLANMODE, patch1[0], trackPar->wR, f1[0], trackPar->wR);
	trackPar->onedForward1R = fftw_create_plan_specific(trackPar->wR,
														FFTW_FORWARD, FFTWPLANMODE, patch1[0], 1, f1[0], 1);
	fprintf(stderr, "5 %i\n", trackPar->wA);
	trackPar->onedForward2 = fftw_create_plan_specific(trackPar->wA,
													   FFTW_FORWARD, FFTWPLANMODE, patch2[0], trackPar->wR, f2[0], trackPar->wR);
	trackPar->onedForward2R = fftw_create_plan_specific(trackPar->wR,
														FFTW_FORWARD, FFTWPLANMODE, patch2[0], 1, f2[0], 1);

	fprintf(stderr, "6 %i %i\n", trackPar->wAa * OS, trackPar->wRa * OS);
	aForward = fftw2d_create_plan(trackPar->wAa * OS, trackPar->wRa * OS,
								  FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "7 %i %i\n", trackPar->wAa * LA * OS, trackPar->wRa * LA * OS);
	aForwardL = fftw2d_create_plan(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS,
								   FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "8 %i %i\n", trackPar->wAa * OS, trackPar->wRa * OS);
	aReverseNoPad = fftw2d_create_plan(trackPar->wAa * OS,
									   trackPar->wRa * OS, FFTW_BACKWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "9 %i %i\n", trackPar->wAa * LA * OS, trackPar->wRa * LA * OS);
	aReverseNoPadL = fftw2d_create_plan(trackPar->wAa * LA * OS,
										trackPar->wRa * LA * OS, FFTW_BACKWARD, FFTWPLANMODE); /* ^^^ */

	fprintf(stderr, "10 %i %i\n", trackPar->wAa, trackPar->wRa);
	aForwardIn = fftw2d_create_plan(trackPar->wAa, trackPar->wRa,
									FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
	fprintf(stderr, "11 %i %i \n", trackPar->wAa * LA, trackPar->wRa * LA);
	aForwardInL = fftw2d_create_plan(trackPar->wAa * LA, trackPar->wRa * LA,
									 FFTW_FORWARD, FFTWPLANMODE); /* ^^^ */
}

/*
  Initialize phase corrections. Malloc spaced, precompute some values.
 */
static void initPhaseCorrect(TrackParams *trackPar)
{
	extern double *thetaD;
	extern double *sinThetaD;
	extern double *cosThetaD;
	extern double *overRange;
	extern fftw_complex **intPatch;
	extern fftw_complex **fftIntPatch;
	extern fftw_complex **intPatchOver;
	extern fftw_complex **fftIntPatchOver;
	double H, Re, lat, range;
	double c1, c2, rc, thetaC, theta;
	int32_t patchSize, osF, nSamps, overSample;
	;
	int32_t i, j;
	osF = trackPar->osF;
	nSamps = trackPar->imageP1.nSlpR;
	if (trackPar->noComplex == TRUE)
		return;
	fprintf(stderr, "init phase Corrections \n");
	if (trackPar->intFlag == TRUE)
	{
		patchSize = max(trackPar->wA / trackPar->intDat.nal, trackPar->wR / trackPar->intDat.nrl);
		patchSize = max(patchSize, INTPATCHSIZE);
		trackPar->intDat.patchSize = patchSize;
		fprintf(stderr, "patchSize %i\n", patchSize);
		fftIntPatch = mallocfftw_complexMat(patchSize, patchSize);
		intPatch = mallocfftw_complexMat(patchSize, patchSize);
		fftIntPatchOver = mallocfftw_complexMat(
			patchSize * trackPar->intDat.nal * osF, patchSize * trackPar->intDat.nrl * osF);
		intPatchOver = mallocfftw_complexMat(
			patchSize * trackPar->intDat.nal * osF, patchSize * trackPar->intDat.nrl * osF);
		for (i = 0; i < patchSize * trackPar->intDat.nal * osF; i++)
			for (j = 0; j < patchSize * trackPar->intDat.nrl * osF; j++)
			{
				fftIntPatchOver[i][j].re = 0;
				fftIntPatchOver[i][j].im = 0;
				intPatchOver[i][j].re = 0;
				intPatchOver[i][j].im = 0;
			}
		trackPar->intDat.forward = fftw2d_create_plan(patchSize,
													  patchSize, FFTW_FORWARD, FFTWPLANMODE);
		trackPar->intDat.backward =
			fftw2d_create_plan(patchSize * trackPar->intDat.nal * osF,
							   trackPar->intDat.nrl * patchSize * osF, FFTW_BACKWARD, FFTWPLANMODE);
	}
	thetaD = (double *)malloc(sizeof(double) * nSamps * osF);
	sinThetaD = (double *)malloc(sizeof(double) * nSamps * osF);
	cosThetaD = (double *)malloc(sizeof(double) * nSamps * osF);
	overRange = (double *)malloc(sizeof(double) * nSamps * osF);
	H = trackPar->imageP1.H;
	lat = trackPar->latc * DTOR;
	Re = earthRadius(lat, EMINOR, EMAJOR) * 1000.0;
	c1 = 2.0 * H * Re + pow(H, 2.0);
	c2 = 2.0 * (Re + H);
	rc = trackPar->imageP1.rn + (double)(nSamps / 2) * trackPar->imageP1.slpR;
	thetaC = acos((pow(rc, 2.0) + c1) / (c2 * rc));
	fprintf(stderr, "rn,H,lat,Re,thetaC,range %f %f %f %f %f %f %i %f\n",
			trackPar->imageP1.rn, H, lat * RTOD, Re, thetaC * RTOD, rc, nSamps, trackPar->imageP1.slpR);

	for (i = 0; i < nSamps * osF; i++)
	{
		range = trackPar->imageP1.rn + (double)i * trackPar->imageP1.slpR / osF;
		theta = acos((pow(range, 2.0) + c1) / (c2 * range));
		thetaD[i] = (theta - thetaC);
		sinThetaD[i] = sin(thetaD[i]);
		cosThetaD[i] = cos(thetaD[i]);
		overRange[i] = 1.0 / range;
	}
}

static void mallocSpace(TrackParams *trackPar)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **fftF1os, **fftF2os;
	extern fftw_complex **fftFa1, **fftFa2;
	extern fftw_complex **psNoPad;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **patch1in, **patch2in;
	extern fftw_complex **cNoPad;
	extern float **cFastOverMag, **cNoPadMag;
	extern float **hanning;
	extern fftw_complex **cFast, **cFastOver;
	extern fftw_complex **img1, **img2;									  /* Detected images for amplitude match*/
	extern fftw_complex **img1in, **img2in;								  /* ^^^ Detected images for amplitude match*/
	extern fftw_complex **fftFa1os, **fftFa1Los, **fftFa2os, **fftFa2Los; /*  ^^^ */
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **caNoPad;
	extern fftw_complex **img1L, **img2L;	  /*images for amplitude match*/
	extern fftw_complex **img1Lin, **img2Lin; /* ^^^ images for amplitude match*/
	extern fftw_complex **psAmpNoPadL;
	extern fftw_complex **fftFa1L, **fftFa2L; /* FFt's */
	extern fftw_complex **caNoPadL;
	extern float **caNoPadMagL;
	extern float **caNoPadMag;
	fftw_complex **fftwTmp;
	strackComplex **ers1Tmp;
	int32_t overSample, i, j;
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
		fprintf(stderr, "Floating point32_t input \n");
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
			for (j = 0; j < imageBuf1.nr; j++)
			{
				fftwTmp[i][j].re = 1;
				fftwTmp[i][j].im = 0;
			}
	}
	else
	{
		fprintf(stderr, "Short int32_t input \n");
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
			for (j = 0; j < imageBuf1.nr; j++)
			{
				ers1Tmp[i][j].r = 1;
				ers1Tmp[i][j].i = 0;
			}
	}
	/*
	  Patches
	*/
	if (trackPar->legacyFlag == TRUE)
		overSample = OSA;
	else
		overSample = 1;
	patch1 = mallocfftw_complexMat(trackPar->wA * overSample, trackPar->wR * overSample);
	patch2 = mallocfftw_complexMat(trackPar->wA * overSample, trackPar->wR * overSample);
	patch1in = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	patch2in = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	img1 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			 /* ^^^ */
	img2 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			 /* ^^^ */
	img1L = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS); /* ^^^ */
	img2L = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS); /* ^^^ */
	img1in = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);					 /* ^^^ */
	img2in = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);					 /* ^^^ */
	img1Lin = mallocfftw_complexMat(trackPar->wAa * LA, trackPar->wRa * LA);		 /* ^^^ */
	img2Lin = mallocfftw_complexMat(trackPar->wAa * LA, trackPar->wRa * LA);		 /* ^^^ */
	for (i = 0; i < trackPar->wAa * OS; i++)
		for (j = 0; j < trackPar->wRa * OS; j++)
		{
			img1[i][j].im = 0.;
			img2[i][j].im = 0.;
		}
	for (i = 0; i < trackPar->wAa * LA * OS; i++)
		for (j = 0; j < trackPar->wRa * LA * OS; j++)
		{
			img1L[i][j].im = 0.;
			img2L[i][j].im = 0.;
		}
	/*
	   One d strided for carrier estimation
	*/
	f1 = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	f2 = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	/*  FFTs	*/
	fftF1 = mallocfftw_complexMat(trackPar->wA * OSA, trackPar->wR * OSA);
	fftF2 = mallocfftw_complexMat(trackPar->wA * OSA, trackPar->wR * OSA);
	fftF1os = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	fftF2os = mallocfftw_complexMat(trackPar->wA, trackPar->wR);
	fftFa1 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			   /* ^^^ */
	fftFa2 = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);			   /* ^^^ */
	fftFa1L = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS); /* ^^^ */
	fftFa2L = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS); /* ^^^ */
	/* ^^^ */
	fftFa1os = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);
	fftFa1Los = mallocfftw_complexMat(trackPar->wAa * LA, trackPar->wRa * LA);
	fftFa2os = mallocfftw_complexMat(trackPar->wAa, trackPar->wRa);
	fftFa2Los = mallocfftw_complexMat(trackPar->wAa * LA, trackPar->wRa * LA);
	/* ^^^ */
	psNoPad = mallocfftw_complexMat(trackPar->wA * OSA, trackPar->wR * OSA);
	psFast = mallocfftw_complexMat(NFAST, NFAST);
	psFastOver = mallocfftw_complexMat(NFAST * NOVER, NFAST * NOVER);
	psAmpNoPad = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);
	psAmpNoPadL = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS);
	zeroComplex(psFastOver, NFAST * NOVER, NFAST * NOVER);
	cFast = mallocfftw_complexMat(NFAST, NFAST);
	cFastOver = mallocfftw_complexMat(NFAST * NOVER, NFAST * NOVER);
	cNoPad = mallocfftw_complexMat(trackPar->wA * OSA, trackPar->wR * OSA);
	caNoPad = mallocfftw_complexMat(trackPar->wAa * OS, trackPar->wRa * OS);
	caNoPadL = mallocfftw_complexMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS);
	cNoPadMag = mallocFloatMat(trackPar->wA * OSA, trackPar->wR * OSA);
	caNoPadMag = mallocFloatMat(trackPar->wAa * OS, trackPar->wRa * OS);
	caNoPadMagL = mallocFloatMat(trackPar->wAa * LA * OS, trackPar->wRa * LA * OS);
	cFastOverMag = mallocFloatMat(NFAST * NOVER, NFAST * NOVER);
	zeroComplex(cFastOver, NFAST * NOVER, NFAST * NOVER);
	zeroFloat(cFastOverMag, NFAST * NOVER, NFAST * NOVER);
	/*      Hanning    */
	hanning = mallocFloatMat(trackPar->wA, trackPar->wR);
	/* Malloc outputs.    */
	trackPar->offR = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->offA = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->corr = mallocFloatMat(trackPar->nA, trackPar->nR);
	trackPar->type = mallocByteMat(trackPar->nA, trackPar->nR);
}

static void trackParInits(TrackParams *trackPar)
{
	char *tBuf, vrtBuf[2048];
	trackPar->nFail = 0.0;
	trackPar->nComplex = 0.0;
	trackPar->nAmp = 0;
	trackPar->nAmpL = 0;
	initPhaseCorrect(trackPar);
	/*
	  Initial offsets initialization
	*/
	if (trackPar->polyShift == FALSE)
	{
		if(readBothOffsetsStrackVrt(&(trackPar->initOffsets), trackPar->initOffsetsFile) == FALSE) {
			// Old method
			trackPar->initOffsets.file = appendSuff(trackPar->initOffsetsFile, ".da", 
				malloc(strlen(trackPar->initOffsetsFile) + 4));
			trackPar->initOffsets.rFile = appendSuff(trackPar->initOffsetsFile, ".dr", 
				malloc(strlen(trackPar->initOffsetsFile) + 4));
			readBothOffsetsStrack(&(trackPar->initOffsets));
		}	
	}
}

/* ************************ Utilities ********************************/

static float absCmpx(fftw_complex x, int32_t noSquareRoot)
{
	float x2;
	x2 = x.re * x.re + x.im * x.im;
	if (noSquareRoot == TRUE)
		return (x2);
	return (sqrt(x2));
}

/* Determine if points (r1,a1) and (r2,a2) with windows wR, wA in their respective image bounds */
static int32_t boundsCheck(int32_t r1, int32_t a1, int32_t r2, int32_t a2, int32_t wR, int32_t wA, TrackParams *trackPar)
{
	int32_t in1, in2;
	in1 = inBounds(r1, a1, wR, wA, trackPar->imageP1.nSlpR, trackPar->imageP1.nSlpA);
	in2 = inBounds(r2, a2, wR, wA, trackPar->imageP2.nSlpR, trackPar->imageP2.nSlpA);
	return (in1 && in2);
}

static void computeHanning(TrackParams *trackPar)
{
	extern float **hanning;
	int32_t i, j, i1, i2, j1, j2;
	double wi, wj;
	fprintf(stderr, "Hanning\n");
	for (i = 0; i < trackPar->wA / 2; i++)
	{
		i1 = trackPar->wA / 2 + i;
		i2 = trackPar->wA / 2 - i - 1;
		wi = 0.5 * (1.0 + cos((2.0 * PI) * (double)i / (double)(trackPar->wA - 1)));
		for (j = 0; j < trackPar->wR / 2; j++)
		{
			j1 = trackPar->wA / 2 + j;
			j2 = trackPar->wA / 2 - j - 1;
			wj = 0.5 * (1.0 + cos((2.0 * PI) * (double)j / (double)(trackPar->wR - 1)));
			hanning[i1][j1] = wi * wj;
			hanning[i2][j2] = wi * wj;
			hanning[i2][j1] = wi * wj;
			hanning[i1][j2] = wi * wj;
		}
	}
}

/*
   Interpolate inital shifts.
 */
static void getShifts(double range, double azimuth, Offsets *offsets, double *dr, double *da)
{
	float t, u, p1, p2, p3, p4;
	float **rimage, **aimage;
	float minvalue, zeroOffset;
	double bn, bp, bSq, deltaZ;
	float alonTrack;

	int32_t i, j;
	minvalue = (float)-LARGEINT * 0.9;
	/*
	  If out of input image, return 0
	*/
	range = (range - offsets->rO) / offsets->deltaR;
	azimuth = (azimuth - offsets->aO) / offsets->deltaA;
	/* Out of bounds */
	if (range < 0.0 || azimuth < 0.0 ||
		(int)range >= (offsets->nr) ||
		(int)azimuth >= (offsets->na))
	{
		*dr = -LARGEINT;
		*da = -LARGEINT;
		return;
	}
	rimage = (float **)offsets->dr;
	aimage = (float **)offsets->da;
	j = (int)range;
	i = (int)azimuth;
	/*
	  Handle border pixels
	*/
	if (j == (offsets->nr - 1) ||
		i == (int)(offsets->na - 1))
	{
		*dr = rimage[i][j];
		*da = aimage[i][j];
		if (*dr < minvalue || *da < minvalue)
		{
			*dr = -LARGEINT;
			*da = -LARGEINT;
		}
		return;
	}
	t = (float)(range - (double)j);
	u = (float)(azimuth - (double)i);
	/* Interp r */
	p1 = rimage[i][j];
	p2 = rimage[i][j + 1];
	p3 = rimage[i + 1][j + 1];
	p4 = rimage[i + 1][j];
	if (p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue)
	{
		*dr = max(max(max(p1, p2), p3), p4);
		if (*dr < minvalue)
		{
			*dr = -LARGEINT;
			*da = -LARGEINT;
			return;
		}
	}
	else
	{
		*dr = (float)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	}
	/* Interp a */
	p1 = aimage[i][j];
	p2 = aimage[i][j + 1];
	p3 = aimage[i + 1][j + 1];
	p4 = aimage[i + 1][j];
	if (p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue)
	{
		*da = max(max(max(p1, p2), p3), p4);
		if (*da < minvalue)
		{
			*dr = -LARGEINT;
			*da = -LARGEINT;
			return;
		}
	}
	else
	{
		*da = (float)((1.0 - t) * (1.0 - u) * p1 + t * (1.0 - u) * p2 + t * u * p3 + (1.0 - t) * u * p4);
	}
}

/* Determine if (r1,r1+wR), (a1,a1+wA) are with thing the bounds (0, sR) and (0,sA) */
static int32_t inBounds(int32_t r1, int32_t a1, int32_t wR, int32_t wA, int32_t sR, int32_t sA)
{
	return (r1 > 0 && (r1 + wR) < sR && a1 > 0 && (a1 + wA) < sA);
}

/*
  NN Interpolate value from mask
 */
static int32_t maskValue(TrackParams *trackPar, int32_t r1, int32_t a1)
{
	int32_t ia1, ir1;
	/* not mask, always try match */
	if (trackPar->maskFlag == FALSE)
		return (1);
	/* otherwise get mask value */
	ia1 = (a1 + trackPar->wA / 2 - trackPar->maskDat.a0) / trackPar->maskDat.nal;
	ir1 = (r1 + trackPar->wR / 2 - trackPar->maskDat.r0) / trackPar->maskDat.nrl;
	if (ia1 < 0 || ia1 >= trackPar->maskDat.na || ir1 < 0 || ir1 >= trackPar->maskDat.nr)
	{
		return (0);
	}
	return ((int)(trackPar->maskDat.mask[ia1][ir1]));
}

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

/*
   Zero pad an array with possible rotation dA, dR (the later 2 are for legacy)
*/
static void zeroPadFFT(fftw_complex **fftPad, int32_t pSA, int32_t pSR, fftw_complex **fftOrig, int32_t oSA, int32_t oSR, int32_t dA, int32_t dR)
{
	int32_t i, j, i1, j1, i2, j2, j2a, ia, ja;
	zeroComplex(fftPad, pSA, pSR);
	/* add non zero elements */
	for (i = 0; i < oSA / 2; i++)
	{
		i1 = pSA - oSA / 2 + i;
		i2 = (oSA / 2 + i + dA) % oSA;
		ia = (i + dA) % oSA;
		j1 = pSR - oSR / 2;
		j2 = oSR / 2 + dR;
		for (j = 0; j < oSR / 2; j++)
		{
			ja = (j + dR) % oSR;
			j2a = j2 % oSR;
			fftPad[i][j] = fftOrig[ia][ja];
			fftPad[i1][j] = fftOrig[i2][ja];
			fftPad[i1][j1] = fftOrig[i2][j2a];
			fftPad[i][j1] = fftOrig[ia][j2a];
			j1++;
			j2++;
		} /* End for j */
	}	  /* End for i */
}

static void zeroComplex(fftw_complex **fft, int32_t na, int32_t nr)
{
	fftw_complex zero;
	int32_t i, j;
	zero.re = 0.;
	zero.im = 0.;
	for (i = 0; i < na; i++)
		for (j = 0; j < nr; j++)
			fft[i][j] = zero;
}

static void zeroFloat(float **x, int32_t na, int32_t nr)
{
	int32_t i, j;
	for (i = 0; i < na; i++)
		for (j = 0; j < nr; j++)
			x[i][j] = 0.;
}

/* ********************************** LEGACY CODE ***************************************/

/*
  Compute  spectrum for oversampled data - legacy case
*/
static void cmpPSNoPad(TrackParams *trackPar)
{
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **psNoPad;
	int32_t i, j, i1;
	FILE *fp1;
	/*  Compute forward fft's	*/
	fftwnd_one(trackPar->cForward, patch1[0], fftF1[0]);
	fftwnd_one(trackPar->cForward, patch2[0], fftF2[0]);
	/*  Zero pad fft for over sampling  spectral shift has been moved to input 11/15/00 ...	*/
	i1 = 0;
	for (i = 0; i < trackPar->wA * OSA; i++)
	{
		for (j = 0; j < trackPar->wR * OSA; j++)
		{
			psNoPad[i][j].re = fftF1[i1][j].re * fftF2[i1][j].re + fftF1[i1][j].im * fftF2[i1][j].im;
			psNoPad[i][j].im = fftF1[i1][j].im * fftF2[i1][j].re - fftF1[i1][j].re * fftF2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
}

/*
  Legacy code to remove estimate any range carrier.
 */
static void estRangeCarrier(TrackParams *trackPar)
{
	extern fftw_complex **patch1in, **patch2in;
	double sp, minS;
	fftw_complex cgrad;
	double angle;
	int32_t iMin;
	int32_t i, j;
	FILE *fp1;
	fftw(trackPar->onedForward1R, trackPar->wA, patch1in[0], 1, trackPar->wR, f1[0], 1, trackPar->wR);
	fftw(trackPar->onedForward2R, trackPar->wA, patch2in[0], 1, trackPar->wR, f2[0], 1, trackPar->wR);
	minS = 2.0e30;
	iMin = 0;
	for (j = 0; j < trackPar->wR; j++)
	{
		sp = 0.;
		for (i = 0; i < trackPar->wA; i++)
		{
			sp += absCmpx(f1[i][j], TRUE);
			sp += absCmpx(f2[i][j], TRUE);
		}
		if (sp < minS)
		{
			iMin = j;
			minS = sp;
		}
	}
	trackPar->rangeShift = (iMin + trackPar->wR / 2) % trackPar->wR;
}

static void azComputeAzShift(int32_t j, int32_t *prevShift, TrackParams *trackPar, double *cosShift, double *sinShift)
{
	double tmpx1, tmpy1;
	if (j > 0 && prevShift >= 0)
	{ /* Avoid wrap around errors */
		tmpx1 = 0.25 * cosShift[trackPar->azShift] + 0.75 * cosShift[*prevShift];
		tmpy1 = 0.25 * sinShift[trackPar->azShift] + 0.75 * sinShift[*prevShift];
		trackPar->azShift = (int)(atan2(tmpy1, tmpx1) / PI * (double)trackPar->wA / 2.0);
		if (trackPar->azShift < 0)
			trackPar->azShift += trackPar->wA;
	}
	*prevShift = trackPar->azShift;
}

/*
  Estimate carrier resulting from non-zero Doppler. This
  is needed to get the oversampling right. The same
  carrier is removed from both patches under the assumption
  that the data were processed to the same doppler.

  Modified 11/14/00 changed to look for min rather
  than max to ensure min gets moved to highend freqs.
*/
static void estCarrier(TrackParams *trackPar)
{
	extern fftw_complex **patch1in, **patch2in;
	extern fftw_complex **f1, **f2;
	double sp, minS;
	fftw_complex cgrad;
	double angle;
	int32_t iMin;
	int32_t i, j;
	fftw(trackPar->onedForward1, trackPar->wR, patch1in[0], trackPar->wR, 1, f1[0], trackPar->wR, 1);
	fftw(trackPar->onedForward2, trackPar->wR, patch2in[0], trackPar->wR, 1, f2[0], trackPar->wR, 1);
	minS = 2.0e30;
	iMin = 0;
	for (i = 0; i < trackPar->wA; i++)
	{
		sp = 0.;
		for (j = 0; j < trackPar->wR; j++)
		{
			sp += f1[i][j].re * f1[i][j].re + f1[i][j].im * f1[i][j].im;
			sp += f2[i][j].re * f2[i][j].re + f2[i][j].im * f2[i][j].im;
		}
		if (sp < minS)
		{
			iMin = i;
			minS = sp;
		}
	}
	trackPar->azShift = (iMin + trackPar->wA / 2) % trackPar->wA;
}

/* ************************ EXPERIMENTAL CODE ********************************/

#define NGAUSS 6
xyData *peakXY = NULL;
float *gaussCorr = NULL;
float *sigGauss = NULL;
float **alpha, **covarm;
float iSigSave = 2.50, jSigSave = 2.5;
/*
Experimental peak fitting with gaussian - not well tested. Correlation peaks not quite right.
*/
static void cmpTrackGauss(TrackParams *trackPar, int32_t *iMax, int32_t *jMax, double *cMax)
{
	extern fftw_complex **fftF1, **fftF2;
	extern fftw_complex **cFast, **cFastOver;
	extern fftw_complex **cNoPad;
	extern float **cFastOverMag, **cNoPadMag;
	extern fftw_complex **patch1, **patch2;
	extern xyData *peakXY;
	extern float *gaussCorr;
	extern float *sigGauss;
	extern float **alpha, **covarm;
	extern float iSigSave, jSigSave;
	double hanningCorrection, corr;
	float a[6]; /* Paramters amp, mux, sigx, muy, sigy, 0 element not used */
	float lasti, lastj;
	float alamada, chisq;		 /* Solver params */
	int32_t i, j, i1, j1, ia[6]; /* ia indicates which params to solve for - should be set a[1:5]=1 */
	int32_t iMax1, jMax1, iData;
	/* This block should only be invoked on the first call to declare global variables */
	if (peakXY == NULL)
	{
		fprintf(stderr, "GaussInit: ");
		peakXY = (xyData *)malloc((NGAUSS * NGAUSS + 1) * sizeof(xyData));
		gaussCorr = vector(1, NGAUSS * NGAUSS);
		sigGauss = vector(1, NGAUSS * NGAUSS);
		covarm = matrix(1, 5, 1, 5);
		alpha = matrix(1, 5, 1, 5);
	}
	for (i = 1; i <= 5; i++)
	{
		ia[i] = 1;
	}
	/* 1.5 removes bias introduced by single hanning window, added 10/19/18 */
	if (trackPar->hanningFlag == TRUE)
		hanningCorrection = 1.5;
	else
		hanningCorrection = 1.0;
	/*   FFT to get correlation function */
	fftwnd_one(trackPar->cReverseNoPad, psNoPad[0], cNoPad[0]);
	/*  Find correlation peak */
	*cMax = 0.;
	findCPeak(cNoPadMag, cNoPad, trackPar->wR * OSA, trackPar->wA * OSA, iMax, jMax, cMax);
	if (*iMax < -999)
		return;
	/*  Load from around peak to oversample;	*/
	iData = 1;
	for (i = 0; i < NGAUSS; i++)
	{
		i1 = *iMax - NGAUSS / 2 + i;
		j1 = *jMax - NGAUSS / 2;
		for (j = 0; j < NGAUSS; j++)
		{
			peakXY[iData].x = i1;
			peakXY[iData].y = j1;
			gaussCorr[iData] = cNoPadMag[i1][j1];
			sigGauss[iData] = 5;
			iData++;
			j1++;
		}
	}
	alamada = -1.;
	a[1] = *cMax;
	a[2] = (float)*iMax;
	a[3] = iSigSave;
	a[4] = (float)*jMax;
	a[5] = jSigSave;
	lasti = 2e9;
	lastj = 2.0e9;
	for (i = 0; i < 5; i++)
	{
		mrqminMod(peakXY, gaussCorr, sigGauss, iData - 1, a, ia, 5, covarm, alpha, &chisq, &myGauss, &alamada);
		if (fabs(a[2] - lasti) < 0.001 && fabs(a[4] - lastj) < 0.01)
			break; /* With NOVER ~ 0.001 */
		lasti = a[2];
		lastj = a[4];
	}
	*iMax = (int)round(a[2] * NOVER);
	*jMax = (int)round(a[4] * NOVER);
	*cMax = (sqrt(a[1])) / sqrt(trackPar->p1 * trackPar->p2);
	*cMax /= trackPar->osF * sqrt(trackPar->wA * trackPar->wR);
}

/*
 Modified version (for xyData) for mrqof from Numerical Recipes in C
 */
#define NRANSI
void mrqcofMod(xyData x[], float y[], float sig[], int32_t ndata, float a[], int32_t ia[],
			   int32_t ma, float **alpha, float beta[], float *chisq,
			   void (*funcs)(xyData, float[], float *, float[], int))
{
	int32_t i, j, k, l, m, mfit = 0;
	float ymod, wt, sig2i, dy, *dyda;
	dyda = vector(1, ma);
	for (j = 1; j <= ma; j++)
		if (ia[j])
			mfit++;
	for (j = 1; j <= mfit; j++)
	{
		for (k = 1; k <= j; k++)
			alpha[j][k] = 0.0;
		beta[j] = 0.0;
	}
	*chisq = 0.0;
	for (i = 1; i <= ndata; i++)
	{
		(*funcs)(x[i], a, &ymod, dyda, ma);
		sig2i = 1.0 / (sig[i] * sig[i]);
		dy = y[i] - ymod;
		for (j = 0, l = 1; l <= ma; l++)
		{
			if (ia[l])
			{
				wt = dyda[l] * sig2i;
				for (j++, k = 0, m = 1; m <= l; m++)
					if (ia[m])
						alpha[j][++k] += wt * dyda[m];
				beta[j] += dy * wt;
			}
		}
		*chisq += dy * dy * sig2i;
	}
	for (j = 2; j <= mfit; j++)
		for (k = 1; k < j; k++)
			alpha[k][j] = alpha[j][k];
	free_vector(dyda, 1, ma);
}

/*
 Modified version (for xyData) of mrqmin from Numerical Recipes in C
 */
void mrqminMod(xyData x[], float y[], float sig[], int32_t ndata, float a[], int32_t ia[],
			   int32_t ma, float **covar, float **alpha, float *chisq,
			   void (*funcs)(xyData, float[], float *, float[], int), float *alamda)
{
	int32_t j, k, l, i;
	static int32_t mfit;
	static float ochisq, *atry, *beta, *da, **oneda;
	if (*alamda < 0.0)
	{
		atry = vector(1, ma);
		beta = vector(1, ma);
		da = vector(1, ma);
		for (mfit = 0, j = 1; j <= ma; j++)
			if (ia[j])
				mfit++;
		oneda = matrix(1, mfit, 1, 1);
		*alamda = 0.001;
		mrqcofMod(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);
		ochisq = (*chisq);
		for (j = 1; j <= ma; j++)
			atry[j] = a[j];
	}
	for (j = 1; j <= mfit; j++)
	{
		for (k = 1; k <= mfit; k++)
			covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
		oneda[j][1] = beta[j];
	}
	gaussj(covar, mfit, oneda, 1);
	for (j = 1; j <= mfit; j++)
		da[j] = oneda[j][1];
	if (*alamda == 0.0)
	{
		covsrt(covar, ma, ia, mfit);
		free_matrix(oneda, 1, mfit, 1, 1);
		free_vector(da, 1, ma);
		free_vector(beta, 1, ma);
		free_vector(atry, 1, ma);
		return;
	}
	for (j = 0, l = 1; l <= ma; l++)
		if (ia[l])
			atry[l] = a[l] + da[++j];
	mrqcofMod(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, funcs);
	if (*chisq < ochisq)
	{
		*alamda *= 0.1;
		ochisq = (*chisq);
		for (j = 1; j <= mfit; j++)
		{
			for (k = 1; k <= mfit; k++)
				alpha[j][k] = covar[j][k];
			beta[j] = da[j];
		}
		for (l = 1; l <= ma; l++)
			a[l] = atry[l];
	}
	else
	{
		*alamda *= 10.0;
		*chisq = ochisq;
	}
}
#undef NRANSI
/* Gaussian function in used by optimization code
   G = f(amp, mux, sigx, muy, sigy) and dG/dparam
 */
static void myGauss(xyData x, float a[], float *y, float dyda[], int32_t na)
{
	float argx, argy, ex, facx, facy;
	float argx1, argy1, ex1, da1, da2, da3, da4, da5;
	int32_t i;
	argx = (x.x - a[2]) / a[3];
	argy = (x.y - a[4]) / a[5];
	ex = exp(-pow(argx, 2) - pow(argy, 2));
	*y = a[1] * ex;
	dyda[1] = ex;
	facx = a[1] * ex * 2 * argx;
	dyda[2] = facx / a[3];
	dyda[3] = facx * argx / a[3];
	facy = a[1] * ex * 2 * argy;
	dyda[4] = facy / a[5];
	dyda[5] = facy * argy / a[5];
}
