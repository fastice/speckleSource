#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
#include <sys/types.h>
#include "time.h"
#include <libgen.h>

#define NOVER 10
#define NFAST 16
#define INTPATCHSIZE 8
#define LA 3
#define BAD 0
#define CMATCH 1
#define AMPMATCH 2
#define AMPMATCHLARGE 3
#define NBUFFERLINES 1000
#define OS 2  /* Amplitude oversample factor */
#define OSA 2

static int getPatches(int r1, int a1, int r2, int a2, FILE *fp1,FILE *fp2, TrackParams *trackPar);
static int getAmpPatches(int r1,int a1,int r2, int a2, FILE *fp1,FILE *fp2, TrackParams *trackPar,int large);
static void findImage2Pos(int r1, int a1, TrackParams *trackPar,   int *r2, int *a2);
static void dumpPatches(fftw_complex **p1,fftw_complex **p2,TrackParams *trackPar);
static void dumpIntPatches(strackComplex **p1, strackComplex **p2,TrackParams *trackPar);
static void fftPlans(TrackParams *trackPar);

static void mallocSpace(TrackParams *trackPar);
static char **mallocByteMat(int nA,int nR);
static float **mallocFloatMat(int nx, int ny);
static fftw_complex **mallocfftw_complexMat(int nx,int ny);
static  strackComplex **mallocErs1ComplexMat(int nA,int nR);
static void cmpPSNoPad(TrackParams *trackPar);
static void overSampleC(TrackParams *trackPar);
static void cmpTrackFast(TrackParams *trackPar, int *iMax, int *jMax, double *cMax);
static void estCarrier(TrackParams *trackPar);
static void estRangeCarrier(TrackParams *trackPar);
static void computeHanning(TrackParams *trackPar);

static void ampMatch(TrackParams *trackPar,int *iMax,int *jMax,
		     double *cMax,int large);
static void ampTrackFast(TrackParams *trackPar, int *iMax, int *jMax, double *cMax,double p1, double p2, int large);
static void findCPeak(float **corrF, fftw_complex **corrNoPad, int wR, int wA, int *iMax,int *jMax,double *cMax);


static void initPhaseCorrect(TrackParams *trackPar);
static void phaseCorrect(TrackParams *trackPar,int r1,int a1);
static void phaseCorrect1(TrackParams *trackPar,int r1,int a1);
static void writeDatFile(TrackParams *trackPar);

static void getShifts(double range,double azimuth,
		      Offsets *offsets,double *dr, double *da);
static void readBothOffsetsStrack( Offsets *offsets);
static int maskValue(TrackParams *trackPar,int r1,int a1);
  
/* ************************ GLOBAL VARIABLES ***************************/
strackComplex *cBuf1,*cBuf2;         /* Input buffers */
fftw_complex *cfBuf1,*cfBuf2;         /* Input buffers */
fftw_complex **fftF1,**fftF2;      /* FFt's */
fftw_complex **patch1, **patch2;
fftw_complex **patch1in, **patch2in; /* ... */
fftw_complex  **psNoPad;
fftw_complex **cFast,**cFastOver;
fftw_complex **cNoPad;
fftw_complex **psFast,**psFastOver;
fftw_complex **f1, **f2;
float **c, **cFastOverMag,**cNoPadMag;
float **hanning;  /* Hanning window */
time_t startTime,lastTime;	
/* Amplitude stuff */
fftw_complex **img1, **img2;       /* Detected images for amplitude match*/
fftw_complex **psAmpNoPad;
fftw_complex **fftFa1,**fftFa2;      /* FFt's */
fftw_complex **caNoPad;
float **caNoPadMag;
/* Double Amplitude stuff */
fftw_complex **img1L, **img2L;  /* Detected images for amplitude match*/
fftw_complex **psAmpNoPadL;
fftw_complex **fftFa1L,**fftFa2L;      /* FFt's */
fftw_complex **caNoPadL;
float **caNoPadMagL;
/* Phase corrections */
double *thetaD;
double *sinThetaD;
double *cosThetaD;
double *overRange;
fftw_complex *bnCorrect;
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
fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
fftw_complex **img1Lin, **img2Lin;/* ^^^ images for amplitude match*/
fftwnd_plan aForwardIn;
fftwnd_plan aForwardInL;
fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
fftwnd_plan cForwardIn; /* ... */
fftw_complex **fftF1os,**fftF2os; /* ... */
int pccount=0;
/* ********************************************************************/




void speckleTrack(TrackParams *trackPar)
{
	extern fftw_complex  **patch1,**patch2, **cFastOver;
	extern float  **cFastOverMag,**cNoPadMag;
	extern fftw_complex **psNoPad;
	extern int pccount;
	FILE *fp1,*fp2;
	FILE *fpR, *fpA, *fpC, *fpT; 
	double tmpx1,tmpy1;
	int i,j;
	int r1,a1,r2,a2;
	int r1a,a1a,r2a,a2a;
	double cMax;
	double wRover2exp, wAover2exp;
	double invNOVER;
	float rShift,aShift,rShift1,aShift1;
	int patchFlag;
	int jMax,iMax;
	int nTot, nCorr,nSkip;
	int good,type;
	int haveData;
	double cAvg;
	double *cosShift,*sinShift;
	int prevShift;
	char *tBuf;
	int maskVal;  /* Flag to determine whether to match */
	int nMask;
	extern time_t startTime,lastTime;	
	fprintf(stderr,"\nSPECKLE TRACKING\n");
	fprintf(stderr,"%i \n",sizeof(fftw_real));
	fprintf(stderr,"%i \n",sizeof(fftw_complex));

	/*
	  Malloc space
	*/
	mallocSpace(trackPar);
	fpR = fopen(trackPar->outFileR,"w");
	fpA = fopen(trackPar->outFileA,"w");
	fpC = fopen(trackPar->outFileC,"w");
	fpT = fopen(trackPar->outFileT,"w");
	sinShift=malloc(trackPar->wA*sizeof(double));
	cosShift=malloc(trackPar->wA*sizeof(double));
	for(i=0; i < trackPar->wA; i++) {
		tmpx1=2.0*PI/(double)(trackPar->wA) * (double)i;
		sinShift[i]=sin(tmpx1);
		cosShift[i]=cos(tmpx1);
	}
	/*
	  compute plan for ffts and hanning window
	*/
	fprintf(stderr,"Plans\n");
	fftPlans(trackPar);
	if (trackPar->hanningFlag == TRUE) {
		fprintf(stderr,"Hanning\n");
		computeHanning(trackPar);
	}
	/*
	  Open input files
	*/
#ifdef HPLARGEFILES    
	fp1 = fopen(trackPar->imageFile1,"r");
	if(fp1 == NULL) error("can't open %s\n",trackPar->imageFile1);
	fp2 = fopen(trackPar->imageFile2,"r");
	if(fp2 == NULL) error("can't open %s\n",trackPar->imageFile2);
#else
	fp1 = openInputFile(trackPar->imageFile1);
	fp2 = openInputFile(trackPar->imageFile2);
#endif

	/*
	  Initializations
	*/ 
	wRover2exp= trackPar->wR/2. * NOVER*OSA;
	wAover2exp= trackPar->wA/2. * NOVER*OSA;
	invNOVER = 1.0/(float)NOVER;
	trackPar->nFail=0.0;
	trackPar->nComplex=0.0;
	trackPar->nAmp=0;
	trackPar->nAmpL=0;
	initPhaseCorrect(trackPar);
	/*
	  Initial offsets initialization 
	*/
	if(trackPar->polyShift==FALSE) {
		tBuf=(char *)malloc(strlen(trackPar->initOffsetsFile)+6);
		fprintf(stderr,"%i\n",strlen(trackPar->initOffsetsFile)+6);
		tBuf[0]='\0';
		tBuf= strcat(tBuf,trackPar->initOffsetsFile);   
		tBuf = strcat(tBuf,".da");    
		trackPar->initOffsets.file=tBuf;
		tBuf=(char *)malloc(strlen(trackPar->initOffsetsFile)+6);
		tBuf[0]='\0';
		tBuf= strcat(tBuf,trackPar->initOffsetsFile);   
		tBuf = strcat(tBuf,".dr");    
		trackPar->initOffsets.rFile=tBuf; 
		readBothOffsetsStrack( &(trackPar->initOffsets));
	}
	/* clear buffers before starting */
	for(i=0; i < trackPar->nA; i++) {
		for(j=0; j < trackPar->nR; j++) {
			trackPar->offR[i][j] = (float)(-LARGEINT);
			trackPar->offA[i][j] =(float)(-LARGEINT);
			trackPar->corr[i][j] = 0;
			trackPar->type[i][j] = BAD;
		}
	}
	/*
	  Loop to do matching 
	*/
	cAvg=0.0;
	good=FALSE;
	prevShift=-1;
	nMask=0;
	nTot=0;
	for(i=0; i < trackPar->nA; i++) {
		lastTime=time(NULL);
		a1= trackPar->aStart + i*trackPar->deltaA - trackPar->wA/2;
		nCorr=0; nSkip=0; cAvg=0.0;
		for(j=0; j < trackPar->nR; j++) {
			rShift=0.0; aShift=0.0; haveData=TRUE;
			r1= trackPar->rStart + j*trackPar->deltaR - trackPar->wR/2;
			/* Find position for second image */
			findImage2Pos(r1,a1,trackPar,&r2,&a2);
			rShift+=(r2-r1); aShift += (a2-a1);
			cMax=0.;
			/* returns 1 if no mask available, 0 or 1 with mask */
			maskVal=maskValue(trackPar,r1*trackPar->scaleFactor,a1*trackPar->scaleFactor);
			if(maskVal > 1) maskVal=1; /* This is because I may put larger values in mask for runcull */
			if(maskVal != 0 && maskVal != 1) fprintf(stderr,"maskVal %i\n",maskVal);
			nMask+=maskVal;
			/* Check Bounds */
			type=BAD;
			if(r1 > 0 && (r1+trackPar->wR) < trackPar->imageP1.nSlpR &&
			   a1 > 0 && (a1+trackPar->wA) < trackPar->imageP1.nSlpA &&
			   r2 > 0 && (r2+trackPar->wR) < trackPar->imageP2.nSlpR &&
			   a2 > 0 && (a2+trackPar->wA) < trackPar->imageP2.nSlpA && maskVal > 0 ) {
				/* Read patches from image */
				haveData=getPatches(r1,a1,r2,a2,fp1,fp2,trackPar);
				/*  Estimate any doppler component in azimuth */
				if(j > 0 && prevShift > 0) {  /* Avoid wrap around errors */
					tmpx1=0.25*cosShift[trackPar->azShift] + 0.75*cosShift[prevShift];
					tmpy1=0.25*sinShift[trackPar->azShift] + 0.75*sinShift[prevShift];
					trackPar->azShift=(int)(atan2(tmpy1,tmpx1)/PI * (double)trackPar->wA/2.0);
					if(trackPar->azShift < 0) trackPar->azShift+=trackPar->wA;
				}
				prevShift=trackPar->azShift;
				/* phase corrections */
				good=FALSE;
				if(trackPar->noComplex==FALSE && haveData==TRUE) {
					/* Phase correction              */
					phaseCorrect1(trackPar, r1, a1);
					/* Complex matching  */
					cmpPSNoPad(trackPar); 
					cmpTrackFast(trackPar,&iMax,&jMax,&cMax);
					if(cMax > 0) {cAvg += cMax; nCorr++;}
					if(cMax > .185 && trackPar->noComplex==FALSE ) {
						/* Update shift with correction */
						rShift1 =((float)jMax-wRover2exp) * invNOVER/OSA;
						aShift1 =((float)iMax-wAover2exp) * invNOVER/OSA;                
						if(fabs((double)rShift1) < 0.125*trackPar->wR*OSA &&
						   fabs((double)aShift1) < 0.125*trackPar->wA*OSA ) {
							good=TRUE; trackPar->nComplex++; type=CMATCH;
						}
					}           
				} else {if(haveData==FALSE) nSkip++;}
			} else good=FALSE;
			/*
			  Amplitude matching 
			*/
			if( good==FALSE && haveData==TRUE && maskVal > 0) {  
				/*	      fprintf(stderr,"."); if(j % 50 == 0) fprintf(stderr,"%i\n",j);*/
				/* Match with small window */
				a1a= trackPar->aStart + i*trackPar->deltaA - trackPar->wAa/2;
				r1a= trackPar->rStart + j*trackPar->deltaR - trackPar->wRa/2;
				r2a=r1a + (r2-r1); a2a=a1a + (a2-a1);
				if(r1a > 0 && (r1a+trackPar->wRa) < trackPar->imageP1.nSlpR &&
				   a1a > 0 && (a1a+trackPar->wAa) < trackPar->imageP1.nSlpA &&
				   r2a > 0 && (r2a+trackPar->wRa) < trackPar->imageP2.nSlpR &&
				   a2a > 0 && (a2a+trackPar->wAa) < trackPar->imageP2.nSlpA ) {
					haveData=getAmpPatches(r1a,a1a,r2a,a2a,fp1,fp2,trackPar,FALSE);
					if(haveData==TRUE) {
						ampMatch(trackPar,&iMax,&jMax,&cMax,FALSE); 
						rShift1 =((float)jMax-trackPar->wRa*0.5*OS* NOVER) * invNOVER/OS;
						aShift1 =((float)iMax-trackPar->wAa*0.5*OS* NOVER) * invNOVER/OS;
					} else { cMax=0.0; nSkip++;}
				} else cMax=0.0;
				if(fabs((double)rShift1) < 0.15*trackPar->wRa*OS &&  fabs((double)aShift1) < 0.15*trackPar->wAa*OS && cMax > .07) { 
					good=TRUE; trackPar->nAmp++; type=AMPMATCH;
				} 
				if(good==FALSE && haveData==TRUE) { /* Match with Large window */
					a1a= trackPar->aStart + i*trackPar->deltaA - trackPar->wAa/2 * LA;
					r1a= trackPar->rStart + j*trackPar->deltaR - trackPar->wRa/2 * LA;
					r2a=r1a + (r2-r1); a2a=a1a + (a2-a1);
					if(r1a > 0 && (r1a+ trackPar->wRa*LA) < trackPar->imageP1.nSlpR &&
					   a1a > 0 && (a1a+ trackPar->wAa*LA) < trackPar->imageP1.nSlpA &&
					   r2a > 0 && (r2a+ trackPar->wRa*LA) < trackPar->imageP2.nSlpR &&
					   a2a > 0 && (a2a+ trackPar->wAa*LA) < trackPar->imageP2.nSlpA ) {
						getAmpPatches(r1a,a1a,r2a,a2a,fp1,fp2,trackPar,TRUE);
						ampMatch(trackPar,&iMax,&jMax,&cMax,TRUE); 
						rShift1 =((float)jMax-trackPar->wRa*0.5 * NOVER*LA*OS) *  invNOVER/OS;
						aShift1 =((float)iMax-trackPar->wAa*0.5 * NOVER*LA*OS) * invNOVER/OS;
						if(fabs((double)rShift1) < 0.175*trackPar->wRa*LA*OS &&
						   fabs((double)aShift1) < 0.175*trackPar->wAa*LA*OS &&
						   cMax > .0) { 
							good=TRUE; trackPar->nAmpL++; type=AMPMATCHLARGE;
						}
					}
				}
			} /* End else for AMp match */
			if(good==TRUE && maskVal > 0) {
				rShift -= rShift1; aShift -= aShift1;
			} else { rShift=-LARGEINT; aShift=-LARGEINT; trackPar->nFail++;}

			/* Save values */

			trackPar->offR[i][j] = rShift;
			trackPar->offA[i][j] = aShift;
			if(rShift > -100000) {trackPar->offR[i][j] *=(float)trackPar->scaleFactor;  trackPar->offA[i][j] *=(float)trackPar->scaleFactor; } 
			trackPar->corr[i][j] = cMax;   /*cMax;*/
			trackPar->type[i][j] = type;   /*cMax;*/
		} /* End for j=0... */

		fwriteBS(trackPar->offR[i],sizeof(float),trackPar->nR,fpR,FLOAT32FLAG);  fflush(fpR);
		fwriteBS(trackPar->offA[i],sizeof(float),trackPar->nR,fpA,FLOAT32FLAG);  fflush(fpA);
		fwriteBS(trackPar->corr[i],sizeof(float),trackPar->nR,fpC,FLOAT32FLAG); fflush(fpC);
		fwriteBS(trackPar->type[i],sizeof(char),trackPar->nR,fpT,BYTEFLAG);  fflush(fpT);
		nTot+=(trackPar->nR - nSkip);
		/*		fprintf(stderr,"nSkip %i \n",nSkip);*/
		if(nCorr > 0) { cAvg=cAvg/(double)(nCorr) ;} else cAvg=0.0;
		fprintf(stderr,
			"\r%6i nTot %i, nMatch %8i %4.1f nComplex %8i %4.1f "
			"nAmp %8i %4.1f (%7i-%7i) nFail %7i %4.1f cAvg %4.2f %i -- %5i (s)",a1,
			nTot,trackPar->nComplex + trackPar->nAmp +trackPar->nAmpL,100.*(double)
			(trackPar->nComplex+trackPar->nAmp+trackPar->nAmpL)/(double)nTot,
			trackPar->nComplex,100.*(double)(trackPar->nComplex)/(double)(nTot),
			trackPar->nAmp+trackPar->nAmpL,    
			100.*(double)(trackPar->nAmp+trackPar->nAmpL)/(double)(nTot),
			trackPar->nAmp,trackPar->nAmpL,
			trackPar->nFail,    100.*(double)(trackPar->nFail)    /(double)(nTot),
			cAvg,nMask,	(int) (time(NULL)-lastTime)	);
	} /* End for i=0... */
	fclose(fpR);
	fclose(fpA);
	fclose(fpC);
	fclose(fpT);
	writeDatFile(trackPar);
}





static void writeDatFile(TrackParams *trackPar)
{
	FILE *fp;
	char geo2[512], *tmp;
	fp = fopen(trackPar->outFileD,"w");
	fprintf(fp,"%i %i %i %i %i %i\n",trackPar->rStart*trackPar->scaleFactor,trackPar->aStart*trackPar->scaleFactor,
		trackPar->nR,trackPar->nA,trackPar->deltaR*trackPar->scaleFactor,trackPar->deltaA*trackPar->scaleFactor);
	/* Add geodats for pair */
	geo2[0]='\0';
	tmp=strcat(geo2,dirname(trackPar->imageFile2));
	tmp=strcat(geo2,trackPar->intGeodat);
	fprintf(fp,"%s %s\n",trackPar->intGeodat,geo2);
}


static void ampMatch(TrackParams *trackPar,int *iMax,int *jMax,
		     double *cMax,int large)
{
	extern fftwnd_plan aForward;
	extern fftwnd_plan aForwardL;
	fftwnd_plan *aFor;
	extern fftw_complex **fftFa1,**fftFa2;
	extern fftw_complex **fftFa1L,**fftFa2L;
	fftw_complex **fa1,**fa2;
	extern fftw_complex **img1,**img2;
	extern fftw_complex **img1L,**img2L;
	fftw_complex **im1,**im2;
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **psAmpNoPadL;
	fftw_complex **psANoPad;
	int wR,wA;
	double avgAmp1,avgAmp2;
	double p1,p2;
	int i,j,i1;
	if(large == TRUE) {
		fa1=fftFa1L; fa2=fftFa2L;
		im1=img1L; im2=img2L;
		psANoPad = psAmpNoPadL;
		wA = trackPar->wAa*LA;
		wR = trackPar->wRa*LA;
		aFor = &aForwardL;
		psANoPad = psAmpNoPadL;
	} else {
		fa1=fftFa1; fa2=fftFa2;
		im1=img1; im2=img2;
		psANoPad = psAmpNoPad;
		wA = trackPar->wAa;
		wR = trackPar->wRa;
		aFor = &aForward;
		psANoPad = psAmpNoPad;
	}
	/* 
	   Step 1: Compute amp and subtract mean
	*/
	avgAmp1 = 0.0;    avgAmp2 = 0.0;
	for(i=0; i < wA*OS; i++) {
		for(j=0; j < wR*OS; j++) {
			im1[i][j].re=sqrt(im1[i][j].re);
			im2[i][j].re=sqrt(im2[i][j].re);
			avgAmp1 += im1[i][j].re;
			avgAmp2 += im2[i][j].re;
		}
	}
	/*
	  Step 2: Subtract mean and est avg pwr
	*/
	avgAmp1=avgAmp1/(double)(wA*OS*wR*OS);
	avgAmp2=avgAmp2/(double)(wA*OS*wR*OS);
	p1=0.0; p2=0.0;
	for(i=0; i < wA*OS; i++) {
		for(j=0; j < wR*OS; j++) {
			if(im1[i][j].re > 1.0e-20) im1[i][j].re -= avgAmp1;
			im2[i][j].re -= avgAmp2;
			p1 += im1[i][j].re * im1[i][j].re;
			p2 += im2[i][j].re * im2[i][j].re;
		}
	}
	p1=p1/(double)(wA*OS*wR*OS);
	p2=p2/(wA*OS*wR*OS);
	/*
	  Step 3: FFT and power spectrum
	*/
	fftwnd_one(*aFor,im1[0],fa1[0]);
	fftwnd_one(*aFor,im2[0],fa2[0]);    
	/*
	  Step 4:Zero pad fft for over sampling
	*/
	/* This now gets done in the oversampling (getamppatches)
	   i1=(int) ( (float)trackPar->azShift * 
	   (float)wA / (float)trackPar->wA   );
	*/
	i1=0;
	for(i=0; i < wA*OS; i++ ) {
		/*         i1= i1 % wA;*/
		for(j=0; j < wR*OS; j++ ) {
			psANoPad[i][j].re = fa1[i1][j].re * fa2[i1][j].re + 
                                fa1[i1][j].im * fa2[i1][j].im; 
			psANoPad[i][j].im = fa1[i1][j].im * fa2[i1][j].re -
                                fa1[i1][j].re * fa2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
	/*
	  Step 5 Compute correlation fuction and find peak
	*/
	ampTrackFast(trackPar,iMax, jMax,cMax,p1,p2,large);
}



static void ampTrackFast(TrackParams *trackPar, int *iMax, int *jMax,
                         double *cMax,double p1, double p2,int large)
{
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **psAmpNoPadL;
	fftw_complex **psANoPad;
	extern fftwnd_plan aReverseNoPad;
	extern fftwnd_plan aReverseNoPadL;
	fftwnd_plan *aReverseNP;
	extern fftw_complex **ps1;
	extern fftw_complex **cFast,**cFastOver;
	extern fftw_complex **caNoPad;
	extern fftw_complex **caNoPadL;
	fftw_complex **caNoP;
	extern fftw_complex **patch1,**patch2;
	extern float **caNoPadMag;
	extern float **caNoPadMagL;
	float **caNoPadM;
	double corr;
	FILE *fp1;
	int i,j,i1,j1;
	int iMax1,jMax1;
	int wA,wR;

	if(large == TRUE) {
		psANoPad = psAmpNoPadL;
		wA = trackPar->wAa*LA; wR = trackPar->wRa*LA;
		aReverseNP = &aReverseNoPadL;
		psANoPad = psAmpNoPadL;
		caNoPadM = caNoPadMagL;
		caNoP = caNoPadL;
	} else {
		psANoPad = psAmpNoPad;
		wA = trackPar->wAa;    wR = trackPar->wRa;
		aReverseNP = &aReverseNoPad;
		psANoPad = psAmpNoPad;
		caNoPadM = caNoPadMag;
		caNoP = caNoPad;
	}
	/* 
	   FFT to get correlation function
	*/
	fftwnd_one(*aReverseNP,psANoPad[0],caNoP[0]);

	findCPeak(caNoPadM,caNoP,wR*OS,wA*OS,iMax,jMax,cMax);

	if(*iMax < -999) return;
	/*  load from around peak to oversample;  */
	for(i=0; i < NFAST; i++) {
		i1 = *iMax -NFAST/2 + i;
		j1 = *jMax - NFAST/2;
		for(j=0; j < NFAST; j++) {
			cFast[i][j].re =caNoPadM[i1][j1];
			cFast[i][j].im = 0.0;
			j1++;
		}
	}
	/* Oversample data */
	overSampleC(trackPar); 
	/*   Find max */
	iMax1=0;
	jMax1=0;
	*cMax=0.0;
	for(i=(NFAST*NOVER)/2-NOVER*2; i <  (NFAST*NOVER)/2 + NOVER*2; i++)
		for(j=(NFAST*NOVER)/2-NOVER*2; j <  (NFAST*NOVER)/2 + NOVER*2; j++) {
			corr = (cFastOver[i][j].re * cFastOver[i][j].re +cFastOver[i][j].im * cFastOver[i][j].im);
			if(corr > *cMax) {*cMax=corr; iMax1=i; jMax1=j;}   
		}
	/* Compute correlation */
	*cMax=  sqrt(*cMax)/(pow((double)(wA*OS),2) *  pow((double)(wR*OS),2) * pow((double) NFAST,2));
	*cMax= *cMax  / sqrt(p1 * p2);
	/*    fprintf(stderr,"%f %f %f %f\n",*cMax,sqrt(p1*p2),p1,p2);*/
	*iMax= (*iMax * NOVER) + iMax1-(NOVER*NFAST)/2;
	*jMax= (*jMax * NOVER) + jMax1-(NOVER*NFAST)/2;
}




/*
  Find peak of no pad correlation function
*/
static void findCPeak(float **corrF, fftw_complex **corrNoPad,
		      int wR, int wA, int *iMax,int *jMax,double *cMax)
{
	int i,j,i1,j1;
	/*
	  Unscramble correlation function
	  can probly omit later
	*/
	*cMax=0.0;

	for(i=0; i < wA/2; i++ ) {
		i1=(wA)/2 + i;
		j1=(wR)/2;
		for(j=0; j < wR/2; j++ ) {
			corrF[i][j] = sqrt(corrNoPad[i1][j1].re * corrNoPad[i1][j1].re +
					   corrNoPad[i1][j1].im * corrNoPad[i1][j1].im);
			if(corrF[i][j] >*cMax)
				{*cMax=corrF[i][j];*iMax=i;*jMax=j; }   

			corrF[i1][j1] = sqrt(corrNoPad[i][j].re * corrNoPad[i][j].re +
					     corrNoPad[i][j].im * corrNoPad[i][j].im);
			if(corrF[i1][j1] > *cMax)
				{*cMax=corrF[i1][j1]; *iMax=i1; *jMax=j1;}   

			corrF[i1][j]  = sqrt(corrNoPad[i][j1].re * corrNoPad[i][j1].re +
					     corrNoPad[i][j1].im * corrNoPad[i][j1].im);
			if(corrF[i1][j] > *cMax) 
				{*cMax=corrF[i1][j]; *iMax=i1; *jMax=j;}   

			corrF[i][j1]  = sqrt(corrNoPad[i1][j].re * corrNoPad[i1][j].re +
					     corrNoPad[i1][j].im * corrNoPad[i1][j].im);
			if(corrF[i][j1] > *cMax) 
				{*cMax=corrF[i][j1]; *iMax=i; *jMax=j1;}   

			j1++;
		}
	}
	/*
	  Check that peak is within range
	*/
	if(*iMax < (NFAST/2) || *iMax > (wA-NFAST/2) ||
	   *jMax < (NFAST/2) || *jMax > (wR-NFAST/2) ) {
		*cMax = -1.0; *jMax=-1000.; *iMax=-1000.;
		return;
	}
}


static int maskValue(TrackParams *trackPar,int r1,int a1)
{
	int ia1,ir1;
	/* not mask, always try match */
	if(trackPar->maskFlag==FALSE) return(1);
	/* otherwise get mask value */
	ia1 = (a1+trackPar->wA/2 - trackPar->maskDat.a0)/trackPar->maskDat.nal;
	ir1 = (r1 +trackPar->wR/2 - trackPar->maskDat.r0)/trackPar->maskDat.nrl;
	if(ia1 < 0 || ia1 >= trackPar->maskDat.na || ir1 < 0 || ir1 >= trackPar->maskDat.nr) { 
		/*		fprintf(stderr,"ia1 ir1 %i %i\n",ia1,ir1); */
		return(0);
	}
	return((int)(trackPar->maskDat.mask[ia1][ir1]));
}



static void phaseCorrect1(TrackParams *trackPar,int r1,int a1)
{
	double x;
	double bn,bp,bsq;
	extern double *sinThetaD;
	extern double *cosThetaD;
	extern double *overRange;
	extern int pccount;
	extern fftw_complex **intPatch;
	extern fftw_complex **fftIntPatch;
	extern fftw_complex **intPatchOver;
	extern fftw_complex **fftIntPatchOver;
	extern fftw_complex **intPatch, **fftIntPatch;
	extern fftw_complex **patch1,**patch2;
	double ptod, deltaX, bnPhase;
	double cr,ci,pr,pi;
	double maxf;
	double tmp;
	double p,scale;
	int i,j,j1,j2,i1,i2;
	int patchSize;
	double lambda;

	int ir1,ia1;  

	patchSize=trackPar->intDat.patchSize;
	x = (a1-trackPar->wA/2-(trackPar->imageP1.nSlpA/2)) /
		(double)(trackPar->imageP1.nSlpA);
	deltaX= 1./(double)(trackPar->imageP1.nSlpA*OSA);
	bn = trackPar->Bn + x * trackPar->dBn + x * x * trackPar->dBnQ;
    
	lambda=trackPar->lambda;
	ptod= - (4.0 * PI / lambda);

	/*
	  Baseline correction
	*/
	for(i=0; i < trackPar->wA*OSA; i++) {
		bp = trackPar->Bp + x * trackPar->dBp + x * x * trackPar->dBpQ;
		bsq = bn*bn + bp*bp;
		for(j=0; j < trackPar->wR*OSA; j++) {
			j1=max(r1*OSA-trackPar->wR*OSA/2,0) + j;
			bnPhase =-bn * sinThetaD[j1] - bp * cosThetaD[j1] + (0.5 * bsq*overRange[j1]);
			bnPhase *= ptod ;
			pr = cos(bnPhase);
			pi = sin(bnPhase);
			cr=patch1[i][j].re;
			ci=patch1[i][j].im;
			patch1[i][j].re = cr * pr - ci * pi;
			patch1[i][j].im = cr * pi + ci * pr;
		}
		x += deltaX;
	} 
	/*
	  Use interferogram to improve match
	*/
	if(trackPar->intFlag == TRUE ) {
		ia1 =(a1+trackPar->wA/2)/trackPar->intDat.nal - trackPar->intDat.patchSize/2;
		ir1 = (r1+trackPar->wR/2)/trackPar->intDat.nrl - trackPar->intDat.patchSize/2;
		if(ia1 > 0 && (ia1+trackPar->intDat.patchSize) < trackPar->intDat.na &&
		   ir1 > 0 && (ir1+trackPar->intDat.patchSize) < trackPar->intDat.nr ) {
			/* Read patch of interferogram data*/          
			for(i=0; i < trackPar->intDat.patchSize; i++) {
				for(j=0; j < trackPar->intDat.patchSize; j++) {
					intPatch[i][j].re=trackPar->intDat.intf[ia1+i][ir1+j].re;
					intPatch[i][j].im=trackPar->intDat.intf[ia1+i][ir1+j].im;
				}
			}
			/* FFT PATCH */
			fftwnd_one(trackPar->intDat.forward,intPatch[0],fftIntPatch[0]);   
			/* Find Peak in interferogram */
			maxf=0.;
			for(i=0; i < trackPar->intDat.patchSize; i++) {
				for(j=0; j < trackPar->intDat.patchSize; j++) {
					p = fftIntPatch[i][j].re *fftIntPatch[i][j].re + fftIntPatch[i][j].im *fftIntPatch[i][j].im;
					maxf=max(maxf,p);
				}
			}
			/* Zero out everything but the peak, which will have the fringe ramp */
			for(i=0; i < trackPar->intDat.patchSize; i++) {
				for(j=0; j < trackPar->intDat.patchSize; j++) {
					p = fftIntPatch[i][j].re *fftIntPatch[i][j].re + fftIntPatch[i][j].im *fftIntPatch[i][j].im;
					if( p < 0.99*maxf ) 
						fftIntPatch[i][j].re=0.0; fftIntPatch[i][j].im =0.0;
				}
			}
			/* OVer Sampling via fft and zero pad */
			for(i=0; i < patchSize/2; i++ ) {
				i1=patchSize*trackPar->intDat.nal*OSA - patchSize/2 + i;
				i2=patchSize/2 + i;
				j1=patchSize*trackPar->intDat.nrl*OSA - patchSize/2;
				j2=patchSize/2;
				for(j=0; j < patchSize/2; j++ ) {
					fftIntPatchOver[i][j].re = fftIntPatch[i][j].re; 
					fftIntPatchOver[i][j].im = fftIntPatch[i][j].im;
					fftIntPatchOver[i1][j].re = fftIntPatch[i2][j].re;
					fftIntPatchOver[i1][j].im = fftIntPatch[i2][j].im;
					fftIntPatchOver[i1][j1].re = fftIntPatch[i2][j2].re;
					fftIntPatchOver[i1][j1].im = fftIntPatch[i2][j2].im;
					fftIntPatchOver[i][j1].re = fftIntPatch[i][j2].re;
					fftIntPatchOver[i][j1].im = fftIntPatch[i][j2].im;
					j1++; j2++;
				} /* End for j */
			} /* End for i */
			fftwnd_one(trackPar->intDat.backward,fftIntPatchOver[0], intPatchOver[0]); 

			/* NOw do correction */
			i1=(patchSize*trackPar->intDat.nal*OSA)/2-trackPar->wA*OSA/2;
			for(i=0; i < trackPar->wA*OSA; i++) {
				j1=(patchSize*trackPar->intDat.nrl*OSA)/2-trackPar->wR*OSA/2;
				for(j=0; j < trackPar->wR*OSA; j++) {
					pr = intPatchOver[i1][j1].re;
					pi = -intPatchOver[i1][j1].im;
					if(fabs((double)pr) < 1.0e-6 && fabs((double)pi) < 1.0e-6) {
						pr=1.0; pi=0.0; scale=1.0;
					} else scale=1.0/sqrt(pr*pr + pi*pi);
					cr=patch1[i][j].re;
					ci=patch1[i][j].im;
					patch1[i][j].re = (cr * pr - ci * pi)*scale;
					patch1[i][j].im = (cr * pi + ci * pr)*scale;
					j1++;
				}
				i1++;
			} 
		} /* End if ia1 > 0 && */
	} /* End if(trackPar->intFlag == TRUE  */
	pccount++;
}


/* SEEMS obsolete */
static void phaseCorrect(TrackParams *trackPar,int r1,int a1)
{
	double x;
	double bn,bp,bsq;
	extern double *sinThetaD;
	extern double *cosThetaD;
	extern double *overRange;
	extern fftw_complex **intPatch;
	extern fftw_complex **fftIntPatch;
	extern fftw_complex **intPatchOver;
	extern fftw_complex **fftIntPatchOver;
	extern fftw_complex **intPatch, **fftIntPatch;
	extern fftw_complex **patch1,**patch2;
	fftw_complex cgrad,cgradr,**intf;
	double phaseA,phaseR;
	double ptod, deltaX, bnPhase;
	double cr,ci,pr,pi;
	double maxf;
	double tmp;
	double p,scale;
	int i,j,j1,j2,i1,i2;
	int patchSize;

	int ir1,ia1;
	patchSize=trackPar->intDat.patchSize;
	x = (a1-trackPar->wA/2-(trackPar->imageP1.nSlpA/2)) /(double)(trackPar->imageP1.nSlpA);
	deltaX= 1./(double)(trackPar->imageP1.nSlpA*OSA);
	bn = trackPar->Bn + x * trackPar->dBn + x * x * trackPar->dBnQ;
	ptod= - (4.0 * PI / trackPar->lambda);
	/*
	  Baseline correction
	*/
	for(i=0; i < trackPar->wA*OSA; i++) {
		bp = trackPar->Bp + x * trackPar->dBp + x * x * trackPar->dBpQ;
		bsq = bn*bn + bp*bp;
		for(j=0; j < trackPar->wR*OSA; j++) {
			j1=max(r1*OSA-trackPar->wR*OSA/2,0) + j;
			bnPhase =-bn * sinThetaD[j1] - bp * cosThetaD[j1] + (0.5 * bsq*overRange[j1]);
			bnPhase *= ptod;
			pr = cos(bnPhase);
			pi = sin(bnPhase);
			cr=patch1[i][j].re;
			ci=patch1[i][j].im;
			patch1[i][j].re = cr * pr - ci * pi;
			patch1[i][j].im = cr * pi + ci * pr;
		}
		x += deltaX;
	}

	/*
	  Use interferogram to improve match
	*/
	if(trackPar->intFlag == TRUE ) {
		ia1 =(a1+trackPar->wA/2)/trackPar->intDat.nal - 
			trackPar->intDat.patchSize/2;
		ir1 = (r1+trackPar->wR/2)/trackPar->intDat.nrl -
			trackPar->intDat.patchSize/2;
		intf=trackPar->intDat.intf;  
		if(ia1 > 0 && (ia1+trackPar->intDat.patchSize) < trackPar->intDat.na &&
		   ir1 > 0 && (ir1+trackPar->intDat.patchSize) < trackPar->intDat.nr ) {
			/* Compute Gradients */    
			cgrad.re=0.0; cgrad.im=0.0; cgradr.re=0.0; cgradr.im=0.0;      
			for(i=ia1; i < (ia1+trackPar->intDat.patchSize-1); i++) {
				for(j=0; j < (ir1+trackPar->intDat.patchSize-1); j++) {
					cgrad.re +=  intf[i][j].re * intf[i+1][j].re + intf[i][j].im * intf[i+1][j].im;
					cgrad.im += -intf[i][j].re * intf[i+1][j].im +intf[i][j].im * intf[i+1][j].re;
					cgradr.re +=  intf[i][j].re * intf[i][j+1].re +intf[i][j].im * intf[i][j+1].im;
					cgradr.im += -intf[i][j].re * intf[i][j+1].im + intf[i][j].im * intf[i][j+1].re;
				}
			}
			phaseA=atan2(cgrad.im,cgrad.re)/(trackPar->intDat.nal*OSA);
			phaseR=atan2(cgradr.im,cgradr.re)/(trackPar->intDat.nrl*OSA);
			/* NOw do correction */
			for(i=0; i < trackPar->wA*OSA; i++) {
				for(j=0; j < trackPar->wR*OSA; j++) {
					pr=cos(phaseA*i+phaseR*j);
					pi=sin(phaseA*i+phaseR*j);
					cr=patch1[i][j].re;
					ci=patch1[i][j].im;
					patch1[i][j].re = (cr * pr - ci * pi);
					patch1[i][j].im = (cr * pi + ci * pr);
					j1++;
				}
				i1++;
			} 
		} /* End if ia1 > 0 && */
	} /* End if(trackPar->intFlag == TRUE  */


}


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
	double H, Re,lat;
	double range;
	double c1,c2,rc,thetaC,theta;
	int nSamps;
	int i,j;
	int patchSize;

	nSamps=trackPar->imageP1.nSlpR;
	fprintf(stderr,"init phase Corrections \n");
	if(trackPar->intFlag==TRUE) {
		patchSize=max(trackPar->wA/trackPar->intDat.nal, trackPar->wR/trackPar->intDat.nrl);
		patchSize=max(patchSize,INTPATCHSIZE);
		trackPar->intDat.patchSize=patchSize;
		fprintf(stderr,"patchSize %i\n",patchSize);
		fftIntPatch = mallocfftw_complexMat(patchSize,patchSize);
		intPatch = mallocfftw_complexMat(patchSize,patchSize);

		fftIntPatchOver = mallocfftw_complexMat(
							patchSize* trackPar->intDat.nal*OSA,patchSize*trackPar->intDat.nrl*OSA);
		intPatchOver = mallocfftw_complexMat(
						     patchSize* trackPar->intDat.nal*OSA,patchSize*trackPar->intDat.nrl*OSA);
		for(i=0; i < patchSize* trackPar->intDat.nal*OSA; i++ )
			for(j=0; j < patchSize* trackPar->intDat.nrl*OSA; j++ ) {
				fftIntPatchOver[i][j].re=0;  fftIntPatchOver[i][j].im=0;
				intPatchOver[i][j].re=0;  intPatchOver[i][j].im=0;
			}

		trackPar->intDat.forward = fftw2d_create_plan(patchSize,
							      patchSize, FFTW_FORWARD, FFTW_MEASURE);
		trackPar->intDat.backward = 
			fftw2d_create_plan(patchSize * trackPar->intDat.nal*OSA,
					   trackPar->intDat.nrl*patchSize*OSA, FFTW_BACKWARD, FFTW_MEASURE);
	}
	thetaD = (double *)malloc(sizeof(double)*nSamps*OSA);
	sinThetaD = (double *)malloc(sizeof(double)*nSamps*OSA);
	cosThetaD = (double *)malloc(sizeof(double)*nSamps*OSA);
	overRange = (double *)malloc(sizeof(double)*nSamps*OSA);
	H=trackPar->imageP1.H;
	lat=trackPar->latc *  DTOR;
	Re=earthRadius(lat,EMINOR,EMAJOR) * 1000.0;  
	c1 = 2.0 * H * Re + pow(H,2.0);
	c2 = 2.0 * (Re + H);
	rc=trackPar->imageP1.rn + (double)(nSamps/2) * trackPar->imageP1.slpR;
	thetaC = acos( (pow(rc,2.0) + c1) / (c2*rc) );
	fprintf(stderr,"rn,H,lat,Re,thetaC,range %f %f %f %f %f %f %i %f\n",trackPar->imageP1.rn ,H,lat*RTOD,Re,thetaC*RTOD,rc, nSamps, trackPar->imageP1.slpR);

	for(i=0; i < nSamps*OSA; i++) {
		range=trackPar->imageP1.rn + (double)i * trackPar->imageP1.slpR/OSA;
		theta = acos( (pow(range,2.0) + c1) / (c2*range) );
		thetaD[i] = (theta - thetaC);
		sinThetaD[i] = sin(thetaD[i]);
		cosThetaD[i] = cos(thetaD[i]);
		overRange[i] = 1.0/range;
	}
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
	extern fftw_complex **f1, **f2;
	extern fftw_complex **patch1in, **patch2in;
	double sp,minS;
	fftw_complex cgrad;
	double angle;
	int iMin;
	int i,j;
	fftw(trackPar->onedForward1,trackPar->wR, patch1in[0],trackPar->wR,1, f1[0],trackPar->wR,1);
	fftw(trackPar->onedForward2,trackPar->wR, patch2in[0],trackPar->wR,1,  f2[0],trackPar->wR,1);
	minS=2.0e30;
	iMin=0;
	for(i=0; i < trackPar->wA; i++) {
		sp=0.;
		for(j=0; j < trackPar->wR; j++) { 
			sp += f1[i][j].re * f1[i][j].re + f1[i][j].im * f1[i][j].im;
			sp += f2[i][j].re * f2[i][j].re + f2[i][j].im * f2[i][j].im;
		}
		if(sp < minS) {
			iMin=i; minS=sp;
		}
	}
	trackPar->azShift=(iMin + trackPar->wA/2) % trackPar->wA;
}


static void estRangeCarrier(TrackParams *trackPar)
{
	extern fftw_complex **f1, **f2;
	extern fftw_complex **patch1in, **patch2in;
	double sp,minS;
	fftw_complex cgrad;
	double angle;
	int iMin;
	int i,j;
	FILE *fp1;
	fftw(trackPar->onedForward1R,trackPar->wA, patch1in[0],1,trackPar->wR,  f1[0],1,trackPar->wR);
	fftw(trackPar->onedForward2R,trackPar->wA, patch2in[0],1,trackPar->wR,  f2[0],1,trackPar->wR);
	minS=2.0e30;
	iMin=0;
	for(j=0; j < trackPar->wR; j++) {
		sp=0.;
		for(i=0; i < trackPar->wA; i++) { 
			sp += f1[i][j].re * f1[i][j].re + f1[i][j].im * f1[i][j].im;
			sp += f2[i][j].re * f2[i][j].re + f2[i][j].im * f2[i][j].im;
		}
		if(sp < minS) {
			iMin=j; minS=sp;
		}
	}
	trackPar->rangeShift=(iMin + trackPar->wR/2) % trackPar->wR;
}


/*
  Compute unpadded spectrum
*/
static void cmpPSNoPad(TrackParams *trackPar)
{
	extern fftw_complex **fftF1,**fftF2;
	extern fftw_complex **patch1,**patch2;
	extern fftw_complex **psNoPad;
	int i,j,i1;
	FILE *fp1;
	/*
	  Compute forward fft's
	*/
	fftwnd_one(trackPar->cForward,patch1[0],fftF1[0]);
	fftwnd_one(trackPar->cForward,patch2[0],fftF2[0]);    
	/*
	  Zero pad fft for over sampling
	  spectral shift has been moved to input 11/15/00 ...
	*/
	i1=0;
	for(i=0; i < trackPar->wA*OSA; i++ ) {
		for(j=0; j < trackPar->wR*OSA; j++ ) {
			psNoPad[i][j].re = fftF1[i1][j].re * fftF2[i1][j].re +
				fftF1[i1][j].im * fftF2[i1][j].im; 
			psNoPad[i][j].im = fftF1[i1][j].im * fftF2[i1][j].re -
				fftF1[i1][j].re * fftF2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
}


static void cmpTrackFast(TrackParams *trackPar, int *iMax, int *jMax,double *cMax)
{
	extern fftw_complex **fftF1,**fftF2;
	extern fftw_complex **ps1;
	extern fftw_complex **cFast,**cFastOver;
	extern fftw_complex **cNoPad;
	extern float  **cFastOverMag,**cNoPadMag;
	extern fftw_complex **patch1,**patch2;
	double hanningCorrection;
	double corr;
	int i,j,i1,j1;
	int iMax1,jMax1;
	FILE *fp1,*fp2;
	/*
		Bias correction for hanning (if used) 
		1.5 removes bias introduced by single hanning window, added 10/19/18 
	*/
	if(trackPar->hanningFlag == TRUE) hanningCorrection=1.5; else hanningCorrection=1.0;
	/* 
	   FFT to get correlation function
	*/
	fftwnd_one(trackPar->cReverseNoPad,psNoPad[0],cNoPad[0]);
	/*
	  Unscramble correlation function
	  can probly omit later
	*/
	*cMax=0.;

	findCPeak(cNoPadMag,cNoPad,trackPar->wR*OSA,trackPar->wA*OSA,iMax,jMax,cMax);

	if(*iMax < -999) return;
	/*
	  load from around peak to oversample;
	*/
	for(i=0; i < NFAST; i++) {
		i1 = *iMax -NFAST/2 + i;
		j1 = *jMax - NFAST/2;
		for(j=0; j < NFAST; j++) {
			cFast[i][j].re =cNoPadMag[i1][j1];
			cFast[i][j].im = 0.0;
			j1++;
		}
	}
	/*
	  Oversample data
	*/
	overSampleC(trackPar); 
	/*
	  Find max
	*/
	iMax1=0;
	jMax1=0;
	*cMax=0.0;
	for(i=(NFAST*NOVER)/2-NOVER*2; i <  (NFAST*NOVER)/2 + NOVER*2; i++)
		for(j=(NFAST*NOVER)/2-NOVER*2; j <  (NFAST*NOVER)/2 + NOVER*2; j++) {
			corr = (cFastOver[i][j].re * cFastOver[i][j].re +	cFastOver[i][j].im * cFastOver[i][j].im);
			if(corr > *cMax) 
				{*cMax=corr; iMax1=i; jMax1=j;}   
		}
       /* Added 10/19/18 for hanning bias */
	*cMax= hanningCorrection* sqrt(*cMax)/(pow((double)(trackPar->wA*OSA),2) * pow((double)(trackPar->wR*OSA),2) * pow((double) NFAST,2));

	*cMax= *cMax  / sqrt(trackPar->p1 * trackPar->p2);
	*iMax=  (*iMax * NOVER) + iMax1-(NOVER*NFAST)/2;
	*jMax= (*jMax * NOVER) + jMax1-(NOVER*NFAST)/2;
}



static void overSampleC(TrackParams *trackPar)
{
	extern fftw_complex **psFast,**psFastOver;
	extern fftw_complex **cFast,**cFastOver;
	int i,j;
	int i1,i2,j1,j2;
	fftwnd_one(trackPar->cForwardFast,cFast[0],psFast[0]);
	/*
	  Zero pad fft for over sampling
	*/
	for(i=0; i < NFAST/2; i++ ) {
		i1=NFAST*NOVER - NFAST/2 +i;
		i2=NFAST/2 + i;
		j1=NFAST*NOVER - NFAST/2;
		j2=NFAST/2;
		for(j=0; j < NFAST/2; j++ ) {
			psFastOver[i][j].re = psFast[i][j].re; 
			psFastOver[i][j].im = psFast[i][j].im;
			psFastOver[i1][j].re = psFast[i2][j].re;
			psFastOver[i1][j].im = psFast[i2][j].im;
			psFastOver[i1][j1].re = psFast[i2][j2].re;
			psFastOver[i1][j1].im = psFast[i2][j2].im;
			psFastOver[i][j1].re = psFast[i][j2].re;
			psFastOver[i][j1].im = psFast[i][j2].im;
			j1++; j2++;
		} /* End for j */
	} /* End for i */
	fftwnd_one(trackPar->cReverseFast,psFastOver[0],cFastOver[0]);
}



static void mallocSpace(TrackParams *trackPar)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;     
	extern strackComplex *cBuf1,*cBuf2;      /* Input buffers */
	extern fftw_complex **fftF1,**fftF2;
	extern fftw_complex **fftF1os,**fftF2os;
	extern fftw_complex **fftFa1,**fftFa2;
	extern fftw_complex  **psNoPad;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **patch1in, **patch2in;
	extern fftw_complex **cNoPad;
	extern float **cFastOverMag, **cNoPadMag; 
	extern float **hanning;
	extern fftw_complex **cFast,**cFastOver;
	extern fftw_complex **f1, **f2;
	extern fftw_complex **img1, **img2;/* Detected images for amplitude match*/
	extern fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
	extern fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **caNoPad;
	extern fftw_complex **img1L, **img2L;/*images for amplitude match*/
	extern fftw_complex **img1Lin, **img2Lin;/* ^^^ images for amplitude match*/
	extern fftw_complex **psAmpNoPadL;
	extern fftw_complex **fftFa1L,**fftFa2L;      /* FFt's */
	extern fftw_complex **caNoPadL;
	fftw_complex **fftwTmp;
	strackComplex **ers1Tmp;
	extern float **caNoPadMagL;
	extern float  **caNoPadMag; 
	int i, j;
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
	if(trackPar->floatFlag ==TRUE) {
		fprintf(stderr,"Floating point input \n");
		imageBuf1.buf=(void **)mallocfftw_complexMat(imageBuf1.na,imageBuf1.nr);
		fftwTmp=(fftw_complex **)imageBuf1.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {fftwTmp[i][j].re=1; fftwTmp[i][j].im=0;}

		imageBuf2.buf=(void**)mallocfftw_complexMat(imageBuf2.na,imageBuf2.nr);
		fftwTmp=(fftw_complex **)imageBuf2.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {fftwTmp[i][j].re=1; fftwTmp[i][j].im=0;}

	}	else {
		fprintf(stderr,"Short int input \n");
		imageBuf1.buf=(void *)mallocErs1ComplexMat(imageBuf1.na,imageBuf1.nr);
		ers1Tmp=(strackComplex **)imageBuf1.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {ers1Tmp[i][j].r=1; ers1Tmp[i][j].i=0;}

		imageBuf2.buf=(void*)mallocErs1ComplexMat(imageBuf2.na,imageBuf2.nr);
		ers1Tmp=(strackComplex **)imageBuf2.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {ers1Tmp[i][j].r=1; ers1Tmp[i][j].i=0;}
	}
	/*
	  Patches
	*/
	patch1= mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	patch2= mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	patch1in= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	patch2in= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	img1= mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	img2= mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	img1L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	img2L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	img1in= mallocfftw_complexMat(trackPar->wAa,trackPar->wRa); /* ^^^ */
	img2in= mallocfftw_complexMat(trackPar->wAa,trackPar->wRa); /* ^^^ */
	img1Lin = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA); /* ^^^ */
	img2Lin = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA); /* ^^^ */
	for(i=0; i < trackPar->wAa*OS; i++) 
		for(j=0; j < trackPar->wRa*OS; j++){img1[i][j].im = 0.; img2[i][j].im = 0.;}
	for(i=0; i < trackPar->wAa*LA*OS; i++) 
		for(j=0; j < trackPar->wRa*LA*OS; j++){
			img1L[i][j].im = 0.; img2L[i][j].im = 0.;}
	/* 
	   One d strided for carrier estimation
	*/
	f1= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	f2= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	/*
	  FFTs
	*/
	fftF1 = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	fftF2 = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	fftF1os = mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	fftF2os = mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	fftFa1 = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	fftFa2 = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	fftFa1L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	fftFa2L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	/* ^^^ */
	fftFa1os =  mallocfftw_complexMat(trackPar->wAa,trackPar->wRa);
	fftFa1Los = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA);
	fftFa2os =  mallocfftw_complexMat(trackPar->wAa,trackPar->wRa);
	fftFa2Los = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA);
	/* ^^^ */

	psNoPad = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	psFast = mallocfftw_complexMat(NFAST,NFAST);
	psFastOver = mallocfftw_complexMat(NFAST*NOVER,NFAST*NOVER);
	psAmpNoPad = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS);
	psAmpNoPadL = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS);

	for(i=0; i < NFAST*NOVER; i++) 
		for(j=0; j < NFAST*NOVER; j++) {
			psFastOver[i][j].re=0.0; psFastOver[i][j].im=0.0;
		}
	/*    c = mallocFloatMat(trackPar->wA*NOVER,trackPar->wR*NOVER); */
	cFast = mallocfftw_complexMat(NFAST,NFAST);
	cFastOver = mallocfftw_complexMat(NFAST*NOVER,NFAST*NOVER);
	cNoPad = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA); 
	caNoPad = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); 
	caNoPadL = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); 
	cNoPadMag = mallocFloatMat(trackPar->wA*OSA,trackPar->wR*OSA); 
	caNoPadMag = mallocFloatMat(trackPar->wAa*OS,trackPar->wRa*OS); 
	caNoPadMagL = mallocFloatMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); 
	cFastOverMag = mallocFloatMat(NFAST*NOVER,NFAST*NOVER);
	for(i=0; i < NFAST*NOVER; i++) 
		for(j=0; j < NFAST*NOVER; j++) {
			cFastOver[i][j].re=0.0; cFastOver[i][j].im=0.0;
			cFastOverMag[i][j]=0.0;
		}
	/*      Hanning    */
	hanning  = mallocFloatMat(trackPar->wA,trackPar->wR);
	/* Input buffers */
	cBuf1=(strackComplex *)
		malloc(sizeof(strackComplex)*max(trackPar->wR,trackPar->wRa*LA));
	cBuf2=(strackComplex *)
		malloc(sizeof(strackComplex)*max(trackPar->wR,trackPar->wRa*LA));

	cfBuf1=(fftw_complex *)
		malloc(sizeof(fftw_complex)*max(trackPar->wR,trackPar->wRa*LA));
	cfBuf2=(fftw_complex *)
		malloc(sizeof(fftw_complex)*max(trackPar->wR,trackPar->wRa*LA));
	/* Malloc outputs.    */
	trackPar->offR =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->offA =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->corr =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->type =mallocByteMat(trackPar->nA,trackPar->nR);

}


static float **mallocFloatMat(int nA,int nR)
{
	float *tmp, **tmp1;
	int i;
	tmp = malloc(nR*nA*sizeof(float));
	tmp1 = (float **)malloc(nA*sizeof(float *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}

static char **mallocByteMat(int nA,int nR)
{
	char *tmp, **tmp1;
	int i;
	tmp = malloc(nR*nA*sizeof(char));
	tmp1 = (char **)malloc(nA*sizeof(char *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}

                         
static fftw_complex **mallocfftw_complexMat(int nA,int nR)
{
	fftw_complex *tmp, **tmp1;
	int i;
	tmp = malloc(nA*nR*sizeof(fftw_complex));
	tmp1 = (fftw_complex **)malloc(nA*sizeof(fftw_complex *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}

static  strackComplex **mallocErs1ComplexMat(int nA,int nR)
{
	strackComplex *tmp, **tmp1;
	int i;
	tmp = malloc(nA*nR*sizeof(strackComplex));
	tmp1 = (strackComplex **)malloc(nA*sizeof(strackComplex *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}


/*
  Read patches
*/
static int getPatches(int r1, int a1, int r2, int a2, FILE *fp1,FILE *fp2, TrackParams *trackPar)
{
#ifdef LARGEFILES
	long long offset1, offset2;
	long long dOff1,dOff2;
#elif defined  HPLARGEFILES
	off_t offset1, offset2;
	off_t dOff1,dOff2;
#else
	long offset1, offset2;
	int dOff1,dOff2;
#endif
	int i,j,i1,j1,i2,j2,ia,ja,j2a;
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;      
	extern strackComplex *cBuf1,*cBuf2; 
	extern fftw_complex *cfBuf1,*cfBuf2; 
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **patch1in, **patch2in;
	extern fftwnd_plan cForwardIn; /* ... */
	double p1,p2;
	strackComplex **ers1Buf1,**ers1Buf2;
	fftw_complex **fftwBuf1,**fftwBuf2;
	fftw_complex zero;
	float scale;
	size_t sSize, wR,s1,s2,a1a,a2a;;
	zero.re=0.0; zero.im=0.0;
	if(a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA) {
		error("getpatches");
		return(FALSE);
	} 
	wR=trackPar->wR;

	if(trackPar->floatFlag ==TRUE) {
		sSize=sizeof(fftw_complex);
		fftwBuf1 = (fftw_complex **)imageBuf1.buf; fftwBuf2 = (fftw_complex **)imageBuf2.buf;
	}else { 
		sSize=sizeof(strackComplex);
		ers1Buf1 = (strackComplex **)imageBuf1.buf; ers1Buf2 = (strackComplex **)imageBuf2.buf;
	}
	/*
	  Load buffers if needed
	*/
	if(a1 < imageBuf1.firstRow || (a1+trackPar->wA) > imageBuf1.lastRow) {
		a1a=max(0,a1-trackPar->wAa*LA);
		s1=trackPar->imageP1.nSlpR*min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES);
		offset1 = ((off_t)a1a * trackPar->imageP1.nSlpR)*(off_t)sSize;
		fseeko(fp1,offset1,SEEK_SET);

		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,FLOAT32FLAG);
		else freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,INT16FLAG);

		imageBuf1.firstRow=a1a;
		imageBuf1.lastRow = a1a + min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES)-1;
	}
	if(a2 < imageBuf2.firstRow || (a2+trackPar->wA) > imageBuf2.lastRow) {
		a2a=max(0,a2-trackPar->wAa*LA);
		s2=trackPar->imageP2.nSlpR*min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES);
		offset2 = ((off_t)a2a * trackPar->imageP2.nSlpR)*(off_t)sSize;
		fseeko(fp2,offset2,SEEK_SET);
		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,FLOAT32FLAG);
		else  freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,INT16FLAG);
		imageBuf2.firstRow=a2a;
		imageBuf2.lastRow = a2a + min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES)-1;
	}

	/*   Loop to input patches    */
	for(i=0; i < trackPar->wA; i++) {
		i1= a1-imageBuf1.firstRow + i;
		i2= a2-imageBuf2.firstRow + i;
		for(j=0; j < trackPar->wR; j++) {
			j1=r1+j;
			j2=r2+j;
			if(trackPar->floatFlag ==TRUE) { 
				patch1in[i][j].re = (fftw_real)fftwBuf1[i1][j1].re;
				patch2in[i][j].re = (fftw_real)fftwBuf2[i2][j2].re;
				patch1in[i][j].im = (fftw_real)fftwBuf1[i1][j1].im;
				patch2in[i][j].im = (fftw_real)fftwBuf2[i2][j2].im;
			} else {
				patch1in[i][j].re = (fftw_real)(ers1Buf1[i1][j1].r);
				patch2in[i][j].re = (fftw_real)(ers1Buf2[i2][j2].r);
				patch1in[i][j].im = (fftw_real)(ers1Buf1[i1][j1].i);
				patch2in[i][j].im = (fftw_real)(ers1Buf2[i2][j2].i);
			}
			if(trackPar->hanningFlag == TRUE) {
				/* removed 10/19/18, compensated correlation with 1.5 scale */
				/*patch1in[i][j].re *= hanning[i][j];*/
				/*patch1in[i][j].im *= hanning[i][j];*/
				patch2in[i][j].re *= hanning[i][j];
				patch2in[i][j].im *= hanning[i][j];
			}
		}
	}
	/* Estimate carrier */
	estCarrier(trackPar);
	estRangeCarrier(trackPar);
	/* Forward FFT */
	fftwnd_one(cForwardIn,patch1in[0],fftF1os[0]);
	fftwnd_one(cForwardIn,patch2in[0],fftF2os[0]);
	/* Zero pad arrays */
	for(i=0; i < trackPar->wA*OSA; i++ ) {
		for(j=0; j < trackPar->wR*OSA; j++ ) {
			fftF1[i][j]=zero;
			fftF2[i][j]=zero;
		}
	}
	for(i=0; i < trackPar->wA/2; i++ ) {
		i1=OSA*trackPar->wA - trackPar->wA/OSA + i;
		i2=(trackPar->wA/OSA + i +  trackPar->azShift) % trackPar->wA;
		ia=(i+ trackPar->azShift) % trackPar->wA;
		j1= trackPar->wR*OSA -  trackPar->wR/OSA;
		j2= trackPar->wR/OSA + trackPar->rangeShift;
		for(j=0; j <  trackPar->wR/2; j++ ) {
			ja=(j+trackPar->rangeShift) % trackPar->wR;
			j2a = j2 % trackPar->wR;
			fftF1[i][j]  = fftF1os[ia][ja]; 
			fftF1[i1][j] = fftF1os[i2][ja];
			fftF1[i1][j1]= fftF1os[i2][j2a];
			fftF1[i][j1] = fftF1os[ia][j2a];
			fftF2[i][j]  = fftF2os[ia][ja]; 
			fftF2[i1][j] = fftF2os[i2][ja];
			fftF2[i1][j1]= fftF2os[i2][j2a];
			fftF2[i][j1] = fftF2os[ia][j2a];
			j1++; j2++;
		} /* End for j */
	} /* End for i */
	/* Invers fft */
	fftwnd_one(trackPar->cReverseNoPad,fftF1[0],patch1[0]);
	fftwnd_one(trackPar->cReverseNoPad,fftF2[0],patch2[0]);
	/* Scale and compute power */
	p1=0.; p2=0.0;
	scale=1.0/(trackPar->wA*trackPar->wR*OSA);
	for(i=0; i < trackPar->wA*OSA; i++ ) {
		for(j=0; j <  trackPar->wR*OSA; j++ ) {
			patch1[i][j].re *= scale;
			patch1[i][j].im *= scale;
			patch2[i][j].re *= scale;
			patch2[i][j].im *= scale; 
			p1 += patch1[i][j].re * patch1[i][j].re +patch1[i][j].im * patch1[i][j].im;
			p2 += patch2[i][j].re * patch2[i][j].re +patch2[i][j].im * patch2[i][j].im;
		}
	}
	if(p1 < 1e-6 || p2 < 1e-6) return(FALSE);

	trackPar->p1 = p1/(double)(trackPar->wA*OSA * trackPar->wR*OSA);
	trackPar->p2 = p2/(double)(trackPar->wA*OSA * trackPar->wR*OSA);
	/* DEBUG   
	   fp1=fopen("patch1","w");    fp2=fopen("patch2","w");
	   fprintf(stderr,"%i\n",sizeof(fftw_complex));
	   fwriteBS((patch1[0]),sizeof(fftw_complex),
	   trackPar->wA*trackPar->wA*OSA*OSA,fp1,FLOAT32FLAG);
	   fwriteBS((patch2[0]),sizeof(fftw_complex),
	   trackPar->wA*trackPar->wA*OSA*OSA,fp2,FLOAT32FLAG);
	   fclose(fp1);fclose(fp2); fp1=fopen("patch1a","w");fp2=fopen("patch2a","w");
	   fprintf(stderr,"%i\n",sizeof(fftw_complex));
	   fwriteBS((patch1in[0]),sizeof(fftw_complex),trackPar->wA*trackPar->wA,fp1,FLOAT32FLAG);
	   fwriteBS((patch2in[0]),sizeof(fftw_complex),trackPar->wA*trackPar->wA,fp2,FLOAT32FLAG);
	   fclose(fp1);fclose(fp2);  exit(-1);
	*/
	return(TRUE);
}




/*
  Read patches
*/
static int getAmpPatches(int r1,int a1,int r2, int a2, FILE *fp1,FILE *fp2,
			 TrackParams *trackPar,int large)
{
#ifdef LARGEFILES
	long long offset1, offset2;
#elif defined  HPLARGEFILES
	off_t offset1, offset2;
#else
	long offset1, offset2;
#endif
	int i,j,i1,j1,i2,j2,ja,j2a; 
	int dOff1,dOff2;
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;     
	extern fftw_complex **img1,**img2;
	extern fftw_complex **img1L,**img2L;
	extern fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
	extern fftw_complex **img1Lin, **img2Lin;/* ^^^ images for amplitude match*/
	extern fftwnd_plan aForwardIn; /* ^^^ */
	extern fftwnd_plan aForwardInL; /* ^^^ */
	extern fftwnd_plan aReverseNoPadL; /* ^^^ */
	extern fftwnd_plan aReverseNoPadL; /* ^^^ */
	extern fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
	fftw_complex **im1,**im2,**im1in,**im2in;
	extern strackComplex *cBuf1,*cBuf2; 
	extern fftw_complex *cfBuf1,*cfBuf2; 
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **fftFa1,**fftFa2;      /* FFt's */
	extern fftw_complex **fftFa1L,**fftFa2L;      /* FFt's */	
	double i1Sum,i2Sum;
	strackComplex **ers1Buf1,**ers1Buf2;
	fftw_complex **fftwBuf1,**fftwBuf2;
	fftw_complex **f1,**f2,**f1a,**f2a;
	int azShift,rangeShift;
	size_t sSize, wRa,wAa,s1,s2,a1a,a2a;
	int ia,i2a;
	fftw_complex zero;
	float scale;
	void *tmp;
	int tmp1;
	FILE *fpJunk;
	zero.re=0.0; zero.im=0.0;
	if(a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA) {
		error("getpatches");
		return(FALSE);
	} 
	if(trackPar->floatFlag ==TRUE) {
		sSize=sizeof(fftw_complex);
		fftwBuf1 = (fftw_complex **)imageBuf1.buf; fftwBuf2 = (fftw_complex **)imageBuf2.buf;
	}else { 
		sSize=sizeof(strackComplex);
		ers1Buf1 = (strackComplex **)imageBuf1.buf; ers1Buf2 = (strackComplex **)imageBuf2.buf;
	}
	if(large == TRUE) {
		wAa=trackPar->wAa*LA;  wRa=trackPar->wRa*LA; im1=img1L; im2=img2L;
		im1in=img1Lin; im2in=img2Lin;
		f1=fftFa1Los;
		f2=fftFa2Los;
		f1a=fftFa1L;
		f2a=fftFa2L;
	} else {              
		wAa=trackPar->wAa;     wRa=trackPar->wRa;    im1=img1; im2=img2; 
		im1in=img1in; im2in=img2in;
		f1=fftFa1os;
		f2=fftFa2os;
		f1a=fftFa1;
		f2a=fftFa2;
	}
	scale=1./(wRa*wAa*wRa*wAa*OS*OS*OS); /* rough guess at scale to avoid fp overflow */
	/*
	  Load buffers if needed
	*/
	if(a1 < imageBuf1.firstRow || (a1+wAa) > imageBuf1.lastRow) {
		/*a1a=max(0,a1-trackPar->wRa*LA);****/
		a1a=max(0,a1-trackPar->wAa*LA);
		s1=trackPar->imageP1.nSlpR*min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES);
		offset1 = ((off_t)a1a * trackPar->imageP1.nSlpR)*(off_t)sSize;
		fseeko(fp1,offset1,SEEK_SET);
		/*	freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,FLOAT32FLAG);*/
		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,FLOAT32FLAG);
		else freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,INT16FLAG);
		imageBuf1.firstRow=a1a;
		imageBuf1.lastRow = a1a + min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES)-1;
	}

	if(a2 < imageBuf2.firstRow || (a2+wAa) > imageBuf2.lastRow) {
		/*a2a=max(0,a2-trackPar->wRa*LA); ****/
		a2a=max(0,a2-trackPar->wAa*LA);
		s2=trackPar->imageP2.nSlpR*min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES);
		offset2 = ((off_t)a2a * trackPar->imageP2.nSlpR)*(off_t)sSize;
		fseeko(fp2,offset2,SEEK_SET);
		/* freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,FLOAT32FLAG);*/
		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,FLOAT32FLAG);
		else  freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,INT16FLAG);
		imageBuf2.firstRow=a2a;
		imageBuf2.lastRow = a2a + min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES)-1;
	}
	/*
	  Get complex data
	*/
	for(i=0; i < wAa; i++) {
		i1= a1-imageBuf1.firstRow + i;
		i2= a2-imageBuf2.firstRow + i;
		for(j=0; j < wRa; j++) {
			j1=r1+j;
			j2=r2+j;
			if(trackPar->floatFlag ==TRUE) {
				im1in[i][j].re = (fftw_real)fftwBuf1[i1][j1].re;
				im2in[i][j].re = (fftw_real)fftwBuf2[i2][j2].re;
				im1in[i][j].im = (fftw_real)fftwBuf1[i1][j1].im;
				im2in[i][j].im = (fftw_real)fftwBuf2[i2][j2].im;
			} else {
				im1in[i][j].re = ((fftw_real)ers1Buf1[i1][j1].r) ;
				im2in[i][j].re = ((fftw_real)ers1Buf2[i2][j2].r) ;
				im1in[i][j].im = ((fftw_real)ers1Buf1[i1][j1].i);
				im2in[i][j].im = ((fftw_real)ers1Buf2[i2][j2].i);
			}
		}
	}
	/* Forward FFT */
	if(large==FALSE) {
		/*fprintf(stderr,"forward fft \n");*/
		fftwnd_one(aForwardIn,im1in[0],f1[0]);
		fftwnd_one(aForwardIn,im2in[0],f2[0]);
	} else {
		fftwnd_one(aForwardInL,im1in[0],f1[0]);
		fftwnd_one(aForwardInL,im2in[0],f2[0]);
	}
	/*
	  Zero pad fft for over sampling
	*/
	for(i=0; i < wAa*OS; i++ ) {
		for(j=0; j < wRa*OS; j++ ) {
			f1a[i][j]=zero;
			f2a[i][j]=zero;
		}
	}
	azShift=(int) ( (float)trackPar->azShift * 	(float)wAa / (float)trackPar->wA   );
	rangeShift=(int) ( (float)trackPar->rangeShift *   (float)wRa / (float)trackPar->wR   );
	for(i=0; i < wAa/2; i++ ) {
		i1=OS*wAa - wAa/OS + i;
		i2=(wAa/OS + i + azShift) % wAa;
		ia=(i+azShift) % wAa;
		j1=wRa*OS - wRa/OS;
		j2= wRa/OS + rangeShift;
		for(j=0; j < wRa/2; j++ ) {
			ja=(j+rangeShift) % wRa;
			j2a = j2 % wRa;
			f1a[i][j]  = f1[ia][ja]; 
			f1a[i1][j] = f1[i2][ja];
			f1a[i1][j1]= f1[i2][j2a];
			f1a[i][j1] = f1[ia][j2a];
			f2a[i][j]  = f2[ia][ja]; 
			f2a[i1][j] = f2[i2][ja];
			f2a[i1][j1]= f2[i2][j2a];
			f2a[i][j1] = f2[ia][j2a];
			j1++; j2++;
		} /* End for j */
	} /* End for i */
	/* Inverse transform */
	if(large==FALSE) {
		fftwnd_one(aReverseNoPad,f1a[0],im1[0]);
		fftwnd_one(aReverseNoPad,f2a[0],im2[0]);
	} else {
		fftwnd_one(aReverseNoPadL,f1a[0],im1[0]);
		fftwnd_one(aReverseNoPadL,f2a[0],im2[0]);
	}
	/*
	  Detect Data 
	*/
	i1Sum=0.0;
	i2Sum=0.0;
	for(i=0; i < wAa*OS; i++) {
		for(j=0; j < wRa*OS; j++) {
			im1[i][j].re=(im1[i][j].re * im1[i][j].re + im1[i][j].im * im1[i][j].im)*scale;
			im2[i][j].re=(im2[i][j].re * im2[i][j].re + im2[i][j].im * im2[i][j].im)*scale;
			i1Sum+=im1[i][j].re;
			i2Sum+=im2[i][j].re;
			im1[i][j].im=0.0;
			im2[i][j].im=0.0;
		}
	}
	/* Return false if either data is just zeros; added 8/11/2016 to avoid zeros in tops images */
	if(i1Sum < 1.0e-6 || i2Sum <1.0e-6) return(FALSE);
	return(TRUE);
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

	trackPar->cForward = fftw2d_create_plan(trackPar->wA*OSA, trackPar->wR*OSA,
						FFTW_FORWARD,FFTW_MEASURE );
	cForwardIn = fftw2d_create_plan(trackPar->wA, trackPar->wR,
					FFTW_FORWARD,FFTW_MEASURE );
	fprintf(stderr,"1 %i %i\n",NFAST,NFAST);

	trackPar->cForwardFast = fftw2d_create_plan(NFAST, NFAST,
						    FFTW_FORWARD,FFTW_MEASURE );
	fprintf(stderr,"2 %i %i\n",trackPar->wA*OSA,trackPar->wR*OSA);
	trackPar->cReverseNoPad = fftw2d_create_plan(trackPar->wA*OSA,
						     trackPar->wR*OSA, FFTW_BACKWARD, FFTW_MEASURE);
	fprintf(stderr,"3 %i %i\n",NOVER * NFAST,NOVER*NFAST );
	trackPar->cReverseFast = fftw2d_create_plan(NOVER * NFAST,
						    NOVER*NFAST, FFTW_BACKWARD, FFTW_MEASURE);
	fprintf(stderr,"4 %i\n",trackPar->wA);
	trackPar->onedForward1 = fftw_create_plan_specific(trackPar->wA,
							   FFTW_FORWARD, FFTW_MEASURE,patch1[0],trackPar->wR,f1[0],trackPar->wR);
	trackPar->onedForward1R = fftw_create_plan_specific(trackPar->wR,
							    FFTW_FORWARD, FFTW_MEASURE,patch1[0],1,f1[0],1);
	fprintf(stderr,"5 %i\n",trackPar->wA);
	trackPar->onedForward2 = fftw_create_plan_specific(trackPar->wA,
							   FFTW_FORWARD, FFTW_MEASURE,patch2[0],trackPar->wR,f2[0],trackPar->wR);
	trackPar->onedForward2R = fftw_create_plan_specific(trackPar->wR,
							    FFTW_FORWARD, FFTW_MEASURE,patch2[0],1,f2[0],1);

	fprintf(stderr,"6 %i %i\n",trackPar->wAa*OS, trackPar->wRa*OS);
	aForward = fftw2d_create_plan(trackPar->wAa*OS, trackPar->wRa*OS,
				      FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */
	fprintf(stderr,"7 %i %i\n", trackPar->wAa*LA*OS, trackPar->wRa*LA*OS);
	aForwardL = fftw2d_create_plan(trackPar->wAa*LA*OS, trackPar->wRa*LA*OS,
				       FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */
	fprintf(stderr,"8 %i %i\n", trackPar->wAa*OS, trackPar->wRa*OS);
	aReverseNoPad = fftw2d_create_plan(trackPar->wAa*OS,
					   trackPar->wRa*OS, FFTW_BACKWARD, FFTW_MEASURE); /* ^^^ */
	fprintf(stderr,"9 %i %i\n",trackPar->wAa*LA*OS, trackPar->wRa*LA*OS);
	aReverseNoPadL = fftw2d_create_plan(trackPar->wAa*LA*OS,
					    trackPar->wRa*LA*OS, FFTW_BACKWARD, FFTW_MEASURE); /* ^^^ */

	fprintf(stderr,"10 %i %i\n",trackPar->wAa, trackPar->wRa);
	aForwardIn = fftw2d_create_plan(trackPar->wAa, trackPar->wRa,
					FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */
	fprintf(stderr,"11 %i %i \n",trackPar->wAa*LA, trackPar->wRa*LA);
	aForwardInL = fftw2d_create_plan(trackPar->wAa*LA, trackPar->wRa*LA,
					 FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */



}


/*
  Find inital shift
*/
static void findImage2Pos(int r1, int a1, TrackParams *trackPar, 
			  int *r2, int *a2)
{
	/*
	  updated 8/21/08 for poly scaling 
	*/
	double aShift,rShift;
	double sc;
	sc=(double)trackPar->scaleFactor;
	if(trackPar->polyShift==TRUE) {
		rShift = trackPar->rShiftPoly[0] + (double)(r1*sc)* trackPar->rShiftPoly[1] +
			(double)(a1*sc)* trackPar->rShiftPoly[2];
		aShift =  trackPar->aShiftPoly[0] + (double)(r1*sc)* trackPar->aShiftPoly[1] +
			(double)(a1*sc)* trackPar->aShiftPoly[2];
	} else {
		getShifts(r1,a1 ,&(trackPar->initOffsets),&rShift,&aShift);
		/* If offset not available use poly - added 6/2/04 */
		if(rShift < (-LARGEINT+1)) {
			rShift = trackPar->rShiftPoly[0] + (double)(r1*sc)* trackPar->rShiftPoly[1] +
				(double)(a1*sc)* trackPar->rShiftPoly[2];
			aShift =  trackPar->aShiftPoly[0] + (double)(r1*sc)* trackPar->aShiftPoly[1] +
				(double)(a1*sc)* trackPar->aShiftPoly[2];
		}
	}
	/* Added 8/21 for tsx scaling */
	rShift /=sc; aShift/=sc; 
	if(rShift >= 0) *r2 = r1 + (int)(rShift + 0.5);
	else *r2 = r1 + (int)(rShift - 0.5);
	if(aShift >= 0) *a2 = a1 + (int)(aShift + 0.5);
	else *a2 = a1 + (int)(aShift - 0.5);
}


static void getShifts(double range,double azimuth,Offsets *offsets,double *dr, double *da)
{
	float t,u;
	float p1,p2,p3,p4;
	float **rimage, **aimage;
	float minvalue;
	float zeroOffset;
	double bn,bp, bSq;
	double deltaZ;
	float alonTrack;
	double imageLength, normAzimuth;
	int i,j;

	minvalue = (float) -LARGEINT * 0.9;
	/*
	  If out of input image, return 0
	*/
	range=(range-offsets->rO)/offsets->deltaR;
	azimuth=(azimuth-offsets->aO)/offsets->deltaA;
	/* Out of bounds */   
	if(range < 0.0 || azimuth < 0.0 ||
	   (int)range >= (offsets->nr) ||
	   (int)azimuth >= (offsets->na) ) {
		*dr=-LARGEINT;  *da=-LARGEINT;
		return;
	}
	rimage = (float **)offsets->dr;  aimage = (float **)offsets->da;
	j = (int)range;    i = (int)azimuth;
	/*
	  Handle border pixels
	*/
	if(j ==  (offsets->nr-1) ||
	   i == (int) (offsets->na-1) ) {
		*dr =  rimage[i][j];  *da =  aimage[i][j];
		if(*dr < minvalue || *da < minvalue) {*dr=-LARGEINT;  *da=-LARGEINT;}
		return;
	}
	t = (float)(range - (double)j);
	u = (float)(azimuth - (double)i);
	/* Interp r */
	p1 = rimage[i][j];        p2 = rimage[i][j+1];
	p3 = rimage[i+1][j+1];    p4 = rimage[i+1][j];
	if(p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue){
		*dr=max(max(max(p1,p2),p3),p4);
		if(*dr < minvalue) {*dr=-LARGEINT; *da=-LARGEINT; return;}
	} else {
		*dr = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +     t * u * p3 +   (1.0 - t) * u* p4);
	}
	/* Interp a */
	p1 = aimage[i][j];        p2 = aimage[i][j+1];
	p3 = aimage[i+1][j+1];    p4 = aimage[i+1][j];
	if(p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue){
		*da=max(max(max(p1,p2),p3),p4);
		if(*da < minvalue) {*dr=-LARGEINT; *da=-LARGEINT; return;}
	} else {
		*da = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +   t * u * p3 +     (1.0 - t) * u* p4);
	}
}


static void computeHanning(TrackParams *trackPar)
{
	extern float **hanning;
	int i,j,i1,i2,j1,j2;
	double wi,wj;
    
	for(i=0; i < trackPar->wA/2; i++ ) {
		i1=trackPar->wA/2+i;
		i2=trackPar->wA/2-i-1;
		wi=0.5*(1.0+cos((2.0*PI)*(double)i/(double)(trackPar->wA -1)));
		for( j=0; j < trackPar->wR/2; j++ ) {
			j1=trackPar->wA/2+j;
			j2=trackPar->wA/2-j-1;
			wj=0.5*(1.0+cos((2.0*PI)*(double)j/(double)(trackPar->wR -1)));
			hanning[i1][j1]=wi*wj;
			hanning[i2][j2]=wi*wj;
			hanning[i2][j1]=wi*wj;
			hanning[i1][j2]=wi*wj;
		}
	}


}




static void readBothOffsetsStrack( Offsets *offsets) 
{
	char *datFile, buf[1024];
	char line[1024];
	int lineCount, eod;
	FILE *fp;
	int rO,aO,nr,na;
	float deltaA,deltaR;
	double c1;
	float *fBuf,*fBufR;
	char *file;
	int i;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=&(buf[0]);
	buf[0]='\0';

	datFile= strcat(datFile,file);   
	datFile = strcat(datFile,".dat");
 
	fp = openInputFile(datFile);
	lineCount=0;
	lineCount=getDataString(fp,lineCount,line,&eod);

	if( sscanf(line,"%i%i%i%i%f%f",&rO,&aO,&nr,&na,&deltaR,&deltaA) != 6)
		error("%s  %i of %s",
		      "readOffsets -- Missing image parameters at line:",
		      lineCount,datFile); 
	fclose(fp);

	offsets->nr=nr;
	offsets->na=na;
	offsets->rO=rO;
	offsets->aO=aO;
	offsets->deltaA = deltaA;
	offsets->deltaR = deltaR;
    

	offsets->da = (float **)malloc(sizeof(float *)*na);
	offsets->dr = (float **)malloc(sizeof(float *)*na);
	/*
	  Malloc buffers
	*/
	fBuf=(float *)malloc(sizeof(float)*na*nr);
	fBufR=(float *)malloc(sizeof(float)*na*nr);
	for(i=0; i < na; i++) {
		offsets->da[i]=&(fBuf[i*nr]); offsets->dr[i]=&(fBufR[i*nr]);
	}
	/*
	  Read files
	*/
	fprintf(stderr,"--- Reading offset file %s\n\n",file);

	fp = openInputFile(offsets->file);
	for(i=0; i < na; i++) freadBS(offsets->da[i],sizeof(float),nr,fp,FLOAT32FLAG);

	fclose(fp);

	fprintf(stderr,"--- Reading offset file %s\n\n",offsets->rFile);

	fp = openInputFile(offsets->rFile);
	for(i=0; i < na; i++) freadBS(offsets->dr[i],sizeof(float),nr,fp,FLOAT32FLAG);
	fclose(fp);
}



/*********************************DEBUG CODE ******************************/



static void dumpPatches(fftw_complex **p1,fftw_complex **p2,TrackParams *trackPar)
{
	FILE *fp1,*fp2;
	fp1=fopen("patch1","w");
	fp2=fopen("patch2","w");
	fprintf(stderr,"%i\n",sizeof(fftw_complex));
	fwriteBS((p1[0]),sizeof(fftw_complex),24*24,fp1,FLOAT32FLAG);
	fwriteBS((p2[0]),sizeof(fftw_complex),24*24,fp2,FLOAT32FLAG);
	exit(-1);
}

static void dumpIntPatches(strackComplex **p1,strackComplex **p2,TrackParams *trackPar)
{
	FILE *fp1,*fp2;
	fp1=fopen("patch1","w");
	fp2=fopen("patch2","w");
	fprintf(stderr,"%i\n",sizeof(fftw_complex));
	fwriteBS((p1[0]),sizeof(strackComplex),24*24,fp1,INT16FLAG);
	fwriteBS((p2[0]),sizeof(strackComplex),24*24,fp2,INT16FLAG);
	exit(-1);
}




/***************************** OBSOLETE CODE *******************************/

/*
  THIS DIDN'T WORK
  Remove carrier resulting from non-zero Doppler. This
  is needed to get the oversampling right. The same
  carrier is removed from both patches under the assumption
  that the data were processed to the same doppler. 

  static void removeCarrier(TrackParams *trackPar)
  {
  extern fftw_complex **f1, **f2;
  extern fftw_complex **patch1, **patch2;
  extern fftw_complex *fMod, *deMod;
  double sp,maxS;
  fftw_complex zero, one;
  int iMax, iPeak;
  int i,j;
  zero.re=0.0;
  zero.im=0.;
  one.re=1.0;
  one.im=0.;

  fftw(trackPar->onedForward1,trackPar->wR,
  patch1[0],trackPar->wR,1,  f1[0],trackPar->wR,1);
  fftw(trackPar->onedForward2,trackPar->wR,
  patch2[0],trackPar->wR,1,  f2[0],trackPar->wR,1);

  maxS=0;
  iMax=0;
  for(i=0; i < trackPar->wA; i++) {
  sp=0.;
  for(j=0; j < trackPar->wR; j++) { 
  sp += f1[i][j].re * f1[i][j].re + f1[i][j].im * f1[i][j].im;
  sp += f2[i][j].re * f2[i][j].re + f2[i][j].im * f2[i][j].im;
  }
  if(sp > maxS) {
  iMax=i; maxS=sp;
  }
  }

  if(iMax > 2 && iMax < (trackPar->wA-2)) {  
  iPeak = trackPar->wA - iMax;
  for(i=0; i < trackPar->wA; i++) 
  if(i == iPeak) fMod[i]=one; else fMod[i]=zero;

  fftw_one(trackPar->deMod,fMod,deMod);


  for(i=0; i < trackPar->wA; i++) {
  for(j=0;  j < trackPar->wR; j++) {
  patch1[i][j].re = patch1[i][j].re * deMod[i].re 
  - patch1[i][j].im * deMod[i].im; 

  patch1[i][j].im = patch1[i][j].re * deMod[i].im 
  + patch1[i][j].im * deMod[i].re; 


  patch2[i][j].re = patch2[i][j].re * deMod[i].re 
  -patch2[i][j].im * deMod[i].im; 

  patch2[i][j].im = patch2[i][j].re * deMod[i].im 
  +patch2[i][j].im * deMod[i].re; 
  }  
  }
  }
  }
*/




/*   static void cmpPS(TrackParams *trackPar)
     {
     extern fftw_complex **fftF1,**fftF2;
     extern fftw_complex **fftR1;
     extern fftw_complex **ps1;
     extern fftw_complex **patch1,**patch2;
     int i,j,i1,i2,j1,j2;


     Zero pad fft for over sampling


     for(i=0; i < trackPar->wA/2; i++ ) {
     i1=trackPar->wA*NOVER - trackPar->wA/2 +i;
     i2=trackPar->wA/2 + i;
     j1=trackPar->wR*NOVER - trackPar->wR/2;
     j2=trackPar->wR/2;
     for(j=0; j < trackPar->wR/2; j++ ) {
     ps1[i][j].re = psNoPad[i][j].re;
     ps1[i1][j].re = psNoPad[i2][j].re;
     ps1[i1][j1].re = psNoPad[i2][j2].re;
     ps1[i][j1].re = psNoPad[i][j2].re;

     ps1[i][j].im = psNoPad[i][j].im;
     ps1[i1][j].im = psNoPad[i2][j].im;
     ps1[i1][j1].im = psNoPad[i2][j2].im;
     ps1[i][j1].im = psNoPad[i][j2].im;
     j1++; j2++;
     } * End for j *
     } * End for i *
     }
*/

/*
  static void cmpTrack(TrackParams *trackPar, int *iMax, int *jMax,
  double *cMax)
  {
  extern fftw_complex **fftF1,**fftF2;
  extern fftw_complex **fftR1;
  extern fftw_complex **ps1;
  extern float **c;
  extern fftw_complex **patch1,**patch2;
  FILE *fp;
  int i,j,i1,i2,j1,j2;
  * 
  FFT to get correlation function
  *
  fftwnd_one(trackPar->cReverse,ps1[0],fftR1[0]);
  *
  Unscramble correlation function
  can probly omit later
  *
  *cMax=0.0;
  for(i=0; i < (trackPar->wA*NOVER)/2; i++ ) {
  i1=(trackPar->wA*NOVER)/2 + i;
  j1=(trackPar->wR*NOVER)/2;
  for(j=0; j < (trackPar->wR*NOVER)/2; j++ ) {

  c[i][j]   = fftR1[i1][j1].re * fftR1[i1][j1].re +
  fftR1[i1][j1].im * fftR1[i1][j1].im;
  if(c[i][j] > *cMax) {*cMax=c[i][j]; *iMax=i; *jMax=j;}   

  c[i1][j1] = fftR1[i][j].re * fftR1[i][j].re +
  fftR1[i][j].im * fftR1[i][j].im;
  if(c[i1][j1] > *cMax) {*cMax=c[i1][j1]; *iMax=i1; *jMax=j1;}   

  c[i1][j]  = fftR1[i][j1].re * fftR1[i][j1].re +
  fftR1[i][j1].im * fftR1[i][j1].im ;
  if(c[i1][j] > *cMax) {*cMax=c[i1][j]; *iMax=i1; *jMax=j;}   

  c[i][j1]  = fftR1[i1][j].re * fftR1[i1][j].re +
  fftR1[i1][j].im * fftR1[i1][j].im ;
  if(c[i][j1] > *cMax) {*cMax=c[i][j1]; *iMax=i; *jMax=j1;}   

  j1++;
  }
  }
  *cMax=
  sqrt(*cMax)/(pow((double)trackPar->wA,2)*pow((double)trackPar->wR,2));
  *cMax=*cMax/sqrt(trackPar->p1 * trackPar->p2);
  }
*/



/*	    	    
		    fprintf(stderr,"%i %i\n",iMax,jMax);
		    fprintf(stderr,"%f %f\n",aShift,rShift);
		    if(i == 1) {
		    dumpPatches(patch1,patch2,trackPar);
		    fp=fopen("junk","w");
		    fwriteBS((psNoPad[0]),sizeof(fftw_complex),trackPar->wR*trackPar->wA,fp,FLOAT32FLAG);
		    fclose(fp);
		    error("stop");
		    }
*/


/* 
   dumpPatches(patch1,patch2,trackPar);
   fp=fopen("junk","w");
   fwriteBS(intPatchOver[0],sizeof(fftw_complex),
   trackPar->intDat.patchSize* trackPar->intDat.patchSize*4*12,fp,FLOAT32FLAG);
   fclose(fp); error("stop");
*/



/* ******************************************************************
   phase gradient estimation
*/
/*    cgrad.re=0;
      cgrad.im=0;
      for(i=1; i < trackPar->wA; i++) {
      for(j=0; j < trackPar->wR; j++) { 
      cgrad.re += f1[i][j].re * f1[i+1][j].re + f1[i][j].im * f1[i+1][j].im;
      cgrad.im += -f1[i][j].re * f1[i+1][j].im + f1[i][j].im * f1[i+1][j].re;
      }
      }
      angle=atan2((double)cgrad.im,(double)cgrad.re)
*/
/* *******************************************************************/


/*      static void phaseCorrect(TrackParams *trackPar,int r1,int a1) */
/*  { */
/*      double x; */
/*      double bn,bp,bsq; */
/*      extern double *sinThetaD; */
/*      extern double *cosThetaD; */
/*      extern double *overRange; */
/*      extern fftw_complex **intPatch; */
/*      extern fftw_complex **fftIntPatch; */
/*      extern fftw_complex **intPatchOver; */
/*      extern fftw_complex **fftIntPatchOver; */
/*      extern fftw_complex **intPatch, **fftIntPatch; */
/*      extern fftw_complex **patch1,**patch2; */
/*      double ptod, deltaX, bnPhase; */
/*      double cr,ci,pr,pi; */
/*      double maxf; */
/*      double tmp; */
/*      double p,scale; */
/*      int i,j,j1,j2,i1,i2; */
/*      int patchSize; */

/*      int ir1,ia1; */
/*      patchSize=trackPar->intDat.patchSize; */
/*      x = (a1-trackPar->wA/2-(trackPar->imageP1.nSlpA/2)) / */
/*           (double)(trackPar->imageP1.nSlpA); */
/*      deltaX= 1./(double)(trackPar->imageP1.nSlpA); */
/*      bn = trackPar->Bn + x * trackPar->dBn + x * x * trackPar->dBnQ; */
/*      ptod= - (4.0 * PI / LAMBDAERS1); */
/*       */
/*          Baseline correction */
/*        */
/*      for(i=0; i < trackPar->wA; i++) { */
/*         bp = trackPar->Bp + x * trackPar->dBp + x * x * trackPar->dBpQ; */
/*         bsq = bn*bn + bp*bp; */
/*         for(j=0; j < trackPar->wR; j++) { */
/*             j1=r1-trackPar->wR/2 + j; */
/*             bnPhase =-bn * sinThetaD[j1] - bp * cosThetaD[j1] +  */
/*             (0.5 * bsq*overRange[j1]); */
/*             bnPhase *= ptod; */
/*             pr = cos(bnPhase); */
/*             pi = sin(bnPhase); */
/*             cr=patch1[i][j].re; */
/*             ci=patch1[i][j].im; */
/*             patch1[i][j].re = cr * pr - ci * pi; */
/*             patch1[i][j].im = cr * pi + ci * pr; */
/*         } */
/*         x += deltaX; */
/*      } */
/*       */
/*        Use interferogram to improve match */
/*       */
/*      if(trackPar->intFlag == TRUE ) { */
/*          ia1 =(a1+trackPar->wA/2)/trackPar->intDat.nal -  */
/*               trackPar->intDat.patchSize/2; */
/*          ir1 = (r1+trackPar->wR/2)/trackPar->intDat.nrl - */
/*               trackPar->intDat.patchSize/2; */
/*          if(ia1 > 0 && (ia1+trackPar->intDat.patchSize) < trackPar->intDat.na && */
/*            ir1 > 0 && (ir1+trackPar->intDat.patchSize) < trackPar->intDat.nr ) { */
/*  	   Read patch            */
/*   	    for(i=0; i < trackPar->intDat.patchSize; i++) { */
/*       	        for(j=0; j < trackPar->intDat.patchSize; j++) { */
/*                     intPatch[i][j].re=trackPar->intDat.intf[ia1+i][ir1+j].re; */
/*                     intPatch[i][j].im=trackPar->intDat.intf[ia1+i][ir1+j].im; */
/*                  } */
/*              } */
/*  	   FFT PATCH  */
/*              fftwnd_one(trackPar->intDat.forward,intPatch[0],fftIntPatch[0]);    */
         
/*             Find Peak */ 
/*              maxf=0.; */
/*              for(i=0; i < trackPar->intDat.patchSize; i++) { */
/*                  for(j=0; j < trackPar->intDat.patchSize; j++) { */
/*                      p = fftIntPatch[i][j].re *fftIntPatch[i][j].re +  */
/*                             fftIntPatch[i][j].im *fftIntPatch[i][j].im; */
/*                      maxf=max(maxf,p); */
/*                  } */
/*              } */
/*               Zero stuff not at peak */
/*              for(i=0; i < trackPar->intDat.patchSize; i++) { */
/*                  for(j=0; j < trackPar->intDat.patchSize; j++) { */
/*                      p = fftIntPatch[i][j].re *fftIntPatch[i][j].re +  */
/*                             fftIntPatch[i][j].im *fftIntPatch[i][j].im; */
/*                      if( p < 0.99*maxf ) { */
/*                          fftIntPatch[i][j].re=0.0; fftIntPatch[i][j].im =0.0; */
/*                      } */
/*                      else fprintf(stderr," %i %i\n",i,j); */
/*                  } */
/*              } */
/*  	    OVer Sampling */
/*              for(i=0; i < patchSize/2; i++ ) { */
/*                  i1=patchSize*trackPar->intDat.nal - patchSize/2 + i; */
/*                  i2=patchSize/2 + i; */
/*                  j1=patchSize*trackPar->intDat.nrl - patchSize/2; */
/*                  j2=patchSize/2; */
/*                  for(j=0; j < patchSize/2; j++ ) { */
/*                      fftIntPatchOver[i][j].re = fftIntPatch[i][j].re;  */
/*                      fftIntPatchOver[i][j].im = fftIntPatch[i][j].im; */
/*                      fftIntPatchOver[i1][j].re = fftIntPatch[i2][j].re; */
/*                      fftIntPatchOver[i1][j].im = fftIntPatch[i2][j].im; */
/*                      fftIntPatchOver[i1][j1].re = fftIntPatch[i2][j2].re; */
/*                      fftIntPatchOver[i1][j1].im = fftIntPatch[i2][j2].im; */
/*                      fftIntPatchOver[i][j1].re = fftIntPatch[i][j2].re; */
/*                      fftIntPatchOver[i][j1].im = fftIntPatch[i][j2].im; */
/*                      j1++; j2++; */
/*                   }  End for j  */
/*               }  End for i  */
/*               fftwnd_one(trackPar->intDat.backward,fftIntPatchOver[0], */
/*                       intPatchOver[0]);  */

/*  	   NOw do correction */
/*               i1=(patchSize*trackPar->intDat.nal)/2-trackPar->wA/2; */
/*               for(i=0; i < trackPar->wA; i++) { */
/*                   j1=(patchSize*trackPar->intDat.nrl)/2-trackPar->wR/2; */
/*                   for(j=0; j < trackPar->wR; j++) { */
/*                      pr = intPatchOver[i1][j1].re; */
/*                      pi = -intPatchOver[i1][j1].im; */
/*                      if(fabs((double)pr) < 1.0e-6 && fabs((double)pi) < 1.0e-6) { */
/*                         pr=1.0; pi=0.0; scale=1.0; */
/*                      } else scale=1.0/sqrt(pr*pr + pi*pi); */
/*                      cr=patch1[i][j].re; */
/*                      ci=patch1[i][j].im; */

/*                      patch1[i][j].re = (cr * pr - ci * pi)*scale; */
/*                      patch1[i][j].im = (cr * pi + ci * pr)*scale; */
		    
/*                      j1++; */
/*                   } */
/*                   i1++; */
/*               }  */
/*  	}  End if ia1 > 0 && */
/*      }  End if(trackPar->intFlag == TRUE   */


/*  } */
