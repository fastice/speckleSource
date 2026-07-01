/*
 * mallocPerThreadArraysW.c
 *
 * Per-thread allocation and deallocation of the strackw work arrays.
 * Called inside a #pragma omp parallel region so each thread initialises
 * its own threadprivate copies.
 *
 * Integration: call before and after the outer i-loop in corrTrackFast():
 *
 *   #pragma omp parallel
 *   { mallocPerThreadArraysW(trackPar); }
 *   for (i = 0; i < trackPar->nA; i++) { ... }
 *   #pragma omp parallel
 *   { freePerThreadArraysW(); }
 */

#include "strackt.h"
#include "clib/standard.h"
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define TP_NOVER 16
#define TP_NFAST 20
#define TP_OS    2

/* Threadprivate globals declared in corrTrackFast.c */
extern fftw_complex **img1, **img2;
extern fftw_complex **img1in, **img2in;
extern fftw_complex **fftFa1, **fftFa2;
extern fftw_complex **fftFa1os, **fftFa2os;
extern fftw_complex **psAmpNoPad;
extern fftw_complex **caNoPad;
extern float **caNoPadMag;
extern fftw_complex **psFast, **psFastOver;
extern fftw_complex **cFast, **cFastOver;
extern double **meanS, **sigmaS, **corrResult;
extern float **dataS, **dataR;
extern double **tmpS1w, **tmpS2w;
extern fftwnd_plan aForward, aReverseNoPad, aForwardIn;
extern fftwnd_plan tp_cForwardFast_w, tp_cReverseFast_w;

/* Must match the threadprivate pragmas in corrTrackFast.c */
#pragma omp threadprivate(img1, img2, img1in, img2in, \
    fftFa1, fftFa2, fftFa1os, fftFa2os, \
    psAmpNoPad, caNoPad, caNoPadMag, \
    psFast, psFastOver, cFast, cFastOver, \
    meanS, sigmaS, corrResult, dataS, dataR, \
    tmpS1w, tmpS2w, \
    aForward, aReverseNoPad, aForwardIn, \
    tp_cForwardFast_w, tp_cReverseFast_w)

/* ------------------------------------------------------------------ */
/* Local helpers: allocate 2-D matrices (contiguous data + row ptrs)  */
/* ------------------------------------------------------------------ */
static fftw_complex **allocCMat(int32_t nA, int32_t nR)
{
    int32_t i;
    fftw_complex *data = malloc((size_t)nA * nR * sizeof(fftw_complex));
    fftw_complex **rows = malloc((size_t)nA * sizeof(fftw_complex *));
    if (!data || !rows)
        error("mallocPerThreadArraysW: out of memory");
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
        error("mallocPerThreadArraysW: out of memory");
    for (i = 0; i < nA; i++)
        rows[i] = &data[i * nR];
    return rows;
}

static double **allocDMat(int32_t nA, int32_t nR)
{
    int32_t i;
    double *data = malloc((size_t)nA * nR * sizeof(double));
    double **rows = malloc((size_t)nA * sizeof(double *));
    if (!data || !rows)
        error("mallocPerThreadArraysW: out of memory");
    for (i = 0; i < nA; i++)
        rows[i] = &data[i * nR];
    return rows;
}

static void freeMC(fftw_complex **m) { if (m) { free(m[0]); free(m); } }
static void freeMF(float **m)        { if (m) { free(m[0]); free(m); } }
static void freeMD(double **m)       { if (m) { free(m[0]); free(m); } }

/* ------------------------------------------------------------------ */
/* Public interface                                                     */
/* ------------------------------------------------------------------ */

void mallocPerThreadArraysW(TrackParams *trackPar)
{
    int32_t wAa = trackPar->wAa;
    int32_t wRa = trackPar->wRa;
    int32_t eA  = trackPar->edgePadA;
    int32_t eR  = trackPar->edgePadR;
    int32_t wA2 = (wAa - 2 * eA) * TP_OS;
    int32_t wR2 = (wRa - 2 * eR) * TP_OS;

    img1   = allocCMat(wAa * TP_OS, wRa * TP_OS);
    img2   = allocCMat(wAa * TP_OS, wRa * TP_OS);
    img1in = allocCMat(wAa, wRa);
    img2in = allocCMat(wAa, wRa);

    fftFa1   = allocCMat(wAa * TP_OS, wRa * TP_OS);
    fftFa2   = allocCMat(wAa * TP_OS, wRa * TP_OS);
    fftFa1os = allocCMat(wAa, wRa);
    fftFa2os = allocCMat(wAa, wRa);

    psAmpNoPad = allocCMat(wAa * TP_OS, wRa * TP_OS);
    caNoPad    = allocCMat(wAa * TP_OS, wRa * TP_OS);
    caNoPadMag = allocFMat(wAa * TP_OS, wRa * TP_OS);

    psFast     = allocCMat(TP_NFAST + 1, TP_NFAST + 1);
    cFast      = allocCMat(TP_NFAST + 1, TP_NFAST + 1);
    cFastOver  = allocCMat((TP_NFAST + 1) * TP_NOVER, (TP_NFAST + 1) * TP_NOVER);

    /* psFastOver is only partially filled in overSampleC (corners); interior
     * must be zeroed so the zero-pad region of the FFT input is correct. */
    psFastOver = allocCMat((TP_NFAST + 1) * TP_NOVER, (TP_NFAST + 1) * TP_NOVER);
    memset(psFastOver[0], 0,
           (size_t)(TP_NFAST + 1) * TP_NOVER * (TP_NFAST + 1) * TP_NOVER * sizeof(fftw_complex));

    meanS      = allocDMat(2 * eA * TP_OS + 1, 2 * eR * TP_OS + 1);
    sigmaS     = allocDMat(2 * eA * TP_OS + 1, 2 * eR * TP_OS + 1);
    corrResult = allocDMat(2 * eA * TP_OS + 1, 2 * eR * TP_OS + 1);
    dataS      = allocFMat(wAa * TP_OS, wRa * TP_OS);
    dataR      = allocFMat(wA2, wR2);
    tmpS1w     = allocDMat(wAa * TP_OS, wRa * TP_OS);
    tmpS2w     = allocDMat(wAa * TP_OS, wRa * TP_OS);

    /* FFTW plan creation is not thread-safe; serialise with a critical section. */
#pragma omp critical(fftw_plan_create)
    {
        aForward      = fftw2d_create_plan(wAa * TP_OS, wRa * TP_OS,
                                           FFTW_FORWARD,  FFTW_ESTIMATE);
        aReverseNoPad = fftw2d_create_plan(wAa * TP_OS, wRa * TP_OS,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);
        aForwardIn    = fftw2d_create_plan(wAa, wRa,
                                           FFTW_FORWARD,  FFTW_ESTIMATE);
        tp_cForwardFast_w = fftw2d_create_plan(TP_NFAST + 1, TP_NFAST + 1,
                                               FFTW_FORWARD,  FFTW_ESTIMATE);
        tp_cReverseFast_w = fftw2d_create_plan(TP_NOVER * (TP_NFAST + 1),
                                               TP_NOVER * (TP_NFAST + 1),
                                               FFTW_BACKWARD, FFTW_ESTIMATE);
    }
}

void freePerThreadArraysW(void)
{
    freeMC(img1);      freeMC(img2);
    freeMC(img1in);    freeMC(img2in);
    freeMC(fftFa1);    freeMC(fftFa2);
    freeMC(fftFa1os);  freeMC(fftFa2os);
    freeMC(psAmpNoPad); freeMC(caNoPad);   freeMF(caNoPadMag);
    freeMC(psFast);    freeMC(psFastOver);
    freeMC(cFast);     freeMC(cFastOver);
    freeMD(meanS);     freeMD(sigmaS);    freeMD(corrResult);
    freeMF(dataS);     freeMF(dataR);
    freeMD(tmpS1w);    freeMD(tmpS2w);
    fftwnd_destroy_plan(aForward);
    fftwnd_destroy_plan(aReverseNoPad);
    fftwnd_destroy_plan(aForwardIn);
    fftwnd_destroy_plan(tp_cForwardFast_w);
    fftwnd_destroy_plan(tp_cReverseFast_w);
}
