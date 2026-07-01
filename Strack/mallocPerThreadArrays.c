/*
 * mallocPerThreadArrays.c
 *
 * Per-thread allocation and deallocation of the speckleTrack work arrays.
 * Called inside a #pragma omp parallel region so each thread initialises its
 * own threadprivate copies.  The fftw_complex** and float** 2-D arrays follow
 * the same row-pointer layout as mallocfftw_complexMat / mallocFloatMat in
 * speckleTrack.c (contiguous data block + separate row-pointer array).
 *
 * Integration: call before and after the outer i-loop in speckleTrack():
 *
 *   #pragma omp parallel
 *   { mallocPerThreadArrays(trackPar); }
 *   for (i = 0; i < trackPar->nA; i++) { ... }
 *   #pragma omp parallel
 *   { freePerThreadArrays(); }
 */

#include "strack.h"
#include <stdlib.h>
#include <omp.h>

/* Constants that match the #define values in speckleTrack.c */
#define TP_NOVER 12
#define TP_NFAST  8
#define TP_OS     2
#define TP_OSA    2
#define TP_LA     2

/* ------------------------------------------------------------------ */
/* Threadprivate globals declared in speckleTrack.c                    */
/* ------------------------------------------------------------------ */
extern fftw_complex **fftF1, **fftF2;
extern fftw_complex **fftF1os, **fftF2os;
extern fftw_complex **fftFa1, **fftFa2;
extern fftw_complex **fftFa1L, **fftFa2L;
extern fftw_complex **fftFa1os, **fftFa1Los, **fftFa2os, **fftFa2Los;
extern fftw_complex **patch1, **patch2;
extern fftw_complex **patch1in, **patch2in;
extern fftw_complex **psNoPad;
extern fftw_complex **psFast, **psFastOver;
extern fftw_complex **cNoPad;
extern fftw_complex **cFast, **cFastOver;
extern fftw_complex **img1, **img2;
extern fftw_complex **img1in, **img2in;
extern fftw_complex **img1L, **img2L;
extern fftw_complex **img1Lin, **img2Lin;
extern fftw_complex **psAmpNoPad;
extern fftw_complex **psAmpNoPadL;
extern fftw_complex **caNoPad;
extern fftw_complex **caNoPadL;
extern fftw_complex **f1, **f2;
extern fftw_complex **intPatch, **fftIntPatch;
extern fftw_complex **intPatchOver, **fftIntPatchOver;
extern float **c, **cFastOverMag, **cNoPadMag;
extern float **caNoPadMag;
extern float **caNoPadMagL;

/* fftwnd plans — file-scope globals and trackPar mirrors, all per-thread */
extern fftwnd_plan aForward, aForwardL;
extern fftwnd_plan aReverseNoPad, aReverseNoPadL;
extern fftwnd_plan aForwardIn, aForwardInL;
extern fftwnd_plan cForwardIn;
extern fftwnd_plan tp_cReverseNoPad, tp_cForwardFast, tp_cReverseFast;
extern fftwnd_plan tp_intForward, tp_intBackward;

/* Must match the threadprivate pragmas in speckleTrack.c */
#pragma omp threadprivate(fftF1, fftF2, fftF1os, fftF2os, \
    patch1, patch2, patch1in, patch2in, \
    psNoPad, cNoPad, cFast, cFastOver, psFast, psFastOver, \
    cFastOverMag, cNoPadMag, c, f1, f2, \
    img1, img2, img1in, img2in, \
    psAmpNoPad, caNoPad, caNoPadMag, \
    img1L, img2L, img1Lin, img2Lin, \
    psAmpNoPadL, fftFa1L, fftFa2L, caNoPadL, caNoPadMagL, \
    fftFa1, fftFa2, \
    fftFa1os, fftFa1Los, fftFa2os, fftFa2Los, \
    intPatch, fftIntPatch, intPatchOver, fftIntPatchOver, \
    aForward, aForwardL, aReverseNoPad, aReverseNoPadL, \
    aForwardIn, aForwardInL, cForwardIn)
#pragma omp threadprivate(tp_cReverseNoPad, tp_cForwardFast, tp_cReverseFast, \
    tp_intForward, tp_intBackward)

/* ------------------------------------------------------------------ */
/* Local helpers: allocate / free 2-D matrices                         */
/* ------------------------------------------------------------------ */
static fftw_complex **allocCMat(int32_t nA, int32_t nR)
{
    int32_t i;
    fftw_complex *data = malloc((size_t)nA * nR * sizeof(fftw_complex));
    fftw_complex **rows = malloc((size_t)nA * sizeof(fftw_complex *));
    if (!data || !rows)
        error("mallocPerThreadArrays: out of memory");
    for (i = 0; i < nA; i++)
        rows[i] = &data[i * nR];
    return rows;
}

static float **allocFMat(int32_t nA, int32_t nR)
{
    int32_t i;
    float *data = malloc((size_t)nA * nR * sizeof(float));
    float **rows = malloc((size_t)nA * sizeof(float *));
    if (!data || !rows)
        error("mallocPerThreadArrays: out of memory");
    for (i = 0; i < nA; i++)
        rows[i] = &data[i * nR];
    return rows;
}

static void freeMC(fftw_complex **m) { if (m) { free(m[0]); free(m); } }
static void freeMF(float **m)        { if (m) { free(m[0]); free(m); } }

/* ------------------------------------------------------------------ */
/* Public interface                                                     */
/* ------------------------------------------------------------------ */

void mallocPerThreadArrays(TrackParams *trackPar)
{
    int32_t wA  = trackPar->wA;
    int32_t wR  = trackPar->wR;
    int32_t wAa = trackPar->wAa;
    int32_t wRa = trackPar->wRa;

    patch1   = allocCMat(wA, wR);
    patch2   = allocCMat(wA, wR);
    patch1in = allocCMat(wA, wR);
    patch2in = allocCMat(wA, wR);

    img1    = allocCMat(wAa * TP_OS,         wRa * TP_OS);
    img2    = allocCMat(wAa * TP_OS,         wRa * TP_OS);
    img1L   = allocCMat(wAa * TP_LA * TP_OS, wRa * TP_LA * TP_OS);
    img2L   = allocCMat(wAa * TP_LA * TP_OS, wRa * TP_LA * TP_OS);
    img1in  = allocCMat(wAa,                 wRa);
    img2in  = allocCMat(wAa,                 wRa);
    img1Lin = allocCMat(wAa * TP_LA,         wRa * TP_LA);
    img2Lin = allocCMat(wAa * TP_LA,         wRa * TP_LA);

    f1 = allocCMat(wA, wR);
    f2 = allocCMat(wA, wR);

    fftF1   = allocCMat(wA * TP_OSA, wR * TP_OSA);
    fftF2   = allocCMat(wA * TP_OSA, wR * TP_OSA);
    fftF1os = allocCMat(wA,          wR);
    fftF2os = allocCMat(wA,          wR);

    fftFa1    = allocCMat(wAa * TP_OS,         wRa * TP_OS);
    fftFa2    = allocCMat(wAa * TP_OS,         wRa * TP_OS);
    fftFa1L   = allocCMat(wAa * TP_LA * TP_OS, wRa * TP_LA * TP_OS);
    fftFa2L   = allocCMat(wAa * TP_LA * TP_OS, wRa * TP_LA * TP_OS);
    fftFa1os  = allocCMat(wAa,                  wRa);
    fftFa2os  = allocCMat(wAa,                  wRa);
    fftFa1Los = allocCMat(wAa * TP_LA,          wRa * TP_LA);
    fftFa2Los = allocCMat(wAa * TP_LA,          wRa * TP_LA);

    psNoPad     = allocCMat(wA * TP_OSA,           wR * TP_OSA);
    psFast      = allocCMat(TP_NFAST,              TP_NFAST);
    psFastOver  = allocCMat(TP_NFAST * TP_NOVER,   TP_NFAST * TP_NOVER);
    cFast       = allocCMat(TP_NFAST,              TP_NFAST);
    cFastOver   = allocCMat(TP_NFAST * TP_NOVER,   TP_NFAST * TP_NOVER);
    psAmpNoPad  = allocCMat(wAa * TP_OS,           wRa * TP_OS);
    psAmpNoPadL = allocCMat(wAa * TP_LA * TP_OS,   wRa * TP_LA * TP_OS);
    cNoPad      = allocCMat(wA * TP_OSA,            wR * TP_OSA);
    caNoPad     = allocCMat(wAa * TP_OS,            wRa * TP_OS);
    caNoPadL    = allocCMat(wAa * TP_LA * TP_OS,    wRa * TP_LA * TP_OS);

    cNoPadMag    = allocFMat(wA * TP_OSA,          wR * TP_OSA);
    caNoPadMag   = allocFMat(wAa * TP_OS,          wRa * TP_OS);
    caNoPadMagL  = allocFMat(wAa * TP_LA * TP_OS,  wRa * TP_LA * TP_OS);
    cFastOverMag = allocFMat(TP_NFAST * TP_NOVER,  TP_NFAST * TP_NOVER);
    c            = allocFMat(1, 1);  /* global declared but unused; placeholder */

    if (trackPar->intFlag == TRUE) {
        int32_t ps  = trackPar->intDat.patchSize;
        int32_t nal = trackPar->intDat.nal;
        int32_t nrl = trackPar->intDat.nrl;
        int32_t osF = trackPar->osF;
        intPatch        = allocCMat(ps, ps);
        fftIntPatch     = allocCMat(ps, ps);
        intPatchOver    = allocCMat(ps * nal * osF, ps * nrl * osF);
        fftIntPatchOver = allocCMat(ps * nal * osF, ps * nrl * osF);
    } else {
        intPatch = fftIntPatch = intPatchOver = fftIntPatchOver = NULL;
    }

    /* Create per-thread fftwnd plans.  FFTW plan creation is not thread-safe
     * (it touches global wisdom state), so serialise with a critical section. */
#pragma omp critical(fftw_plan_create)
    {
        cForwardIn    = fftw2d_create_plan(wA,                 wR,                 FFTW_FORWARD,  FFTW_ESTIMATE);
        aForward      = fftw2d_create_plan(wAa * TP_OS,        wRa * TP_OS,        FFTW_FORWARD,  FFTW_ESTIMATE);
        aForwardL     = fftw2d_create_plan(wAa * TP_LA * TP_OS,wRa * TP_LA * TP_OS,FFTW_FORWARD, FFTW_ESTIMATE);
        aReverseNoPad = fftw2d_create_plan(wAa * TP_OS,        wRa * TP_OS,        FFTW_BACKWARD, FFTW_ESTIMATE);
        aReverseNoPadL= fftw2d_create_plan(wAa * TP_LA * TP_OS,wRa * TP_LA * TP_OS,FFTW_BACKWARD,FFTW_ESTIMATE);
        aForwardIn    = fftw2d_create_plan(wAa,                wRa,                FFTW_FORWARD,  FFTW_ESTIMATE);
        aForwardInL   = fftw2d_create_plan(wAa * TP_LA,        wRa * TP_LA,        FFTW_FORWARD,  FFTW_ESTIMATE);
        tp_cReverseNoPad = fftw2d_create_plan(wA * TP_OSA,     wR * TP_OSA,        FFTW_BACKWARD, FFTW_ESTIMATE);
        tp_cForwardFast  = fftw2d_create_plan(TP_NFAST,        TP_NFAST,           FFTW_FORWARD,  FFTW_ESTIMATE);
        tp_cReverseFast  = fftw2d_create_plan(TP_NOVER*TP_NFAST,TP_NOVER*TP_NFAST, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (trackPar->intFlag == TRUE) {
            int32_t ps  = trackPar->intDat.patchSize;
            int32_t nal = trackPar->intDat.nal;
            int32_t nrl = trackPar->intDat.nrl;
            int32_t osF = trackPar->osF;
            tp_intForward  = fftw2d_create_plan(ps,           ps,           FFTW_FORWARD,  FFTW_ESTIMATE);
            tp_intBackward = fftw2d_create_plan(ps * nal * osF, ps * nrl * osF, FFTW_BACKWARD, FFTW_ESTIMATE);
        } else {
            tp_intForward = tp_intBackward = NULL;
        }
    }
}

void freePerThreadArrays(void)
{
    freeMC(patch1);    freeMC(patch2);
    freeMC(patch1in);  freeMC(patch2in);
    freeMC(img1);      freeMC(img2);
    freeMC(img1L);     freeMC(img2L);
    freeMC(img1in);    freeMC(img2in);
    freeMC(img1Lin);   freeMC(img2Lin);
    freeMC(f1);        freeMC(f2);
    freeMC(fftF1);     freeMC(fftF2);
    freeMC(fftF1os);   freeMC(fftF2os);
    freeMC(fftFa1);    freeMC(fftFa2);
    freeMC(fftFa1L);   freeMC(fftFa2L);
    freeMC(fftFa1os);  freeMC(fftFa2os);
    freeMC(fftFa1Los); freeMC(fftFa2Los);
    freeMC(psNoPad);
    freeMC(psFast);    freeMC(psFastOver);
    freeMC(cFast);     freeMC(cFastOver);
    freeMC(psAmpNoPad);   freeMC(psAmpNoPadL);
    freeMC(cNoPad);       freeMC(caNoPad);    freeMC(caNoPadL);
    freeMF(cNoPadMag);    freeMF(caNoPadMag); freeMF(caNoPadMagL);
    freeMF(cFastOverMag); freeMF(c);
    freeMC(intPatch);     freeMC(fftIntPatch);
    freeMC(intPatchOver); freeMC(fftIntPatchOver);
    fftwnd_destroy_plan(cForwardIn);
    fftwnd_destroy_plan(aForward);      fftwnd_destroy_plan(aForwardL);
    fftwnd_destroy_plan(aReverseNoPad); fftwnd_destroy_plan(aReverseNoPadL);
    fftwnd_destroy_plan(aForwardIn);    fftwnd_destroy_plan(aForwardInL);
    fftwnd_destroy_plan(tp_cReverseNoPad);
    fftwnd_destroy_plan(tp_cForwardFast);
    fftwnd_destroy_plan(tp_cReverseFast);
    if (tp_intForward)  fftwnd_destroy_plan(tp_intForward);
    if (tp_intBackward) fftwnd_destroy_plan(tp_intBackward);
}
