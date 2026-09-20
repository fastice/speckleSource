#ifndef STRACK_FFT_H
#define STRACK_FFT_H
#include <stdint.h>
#include <stddef.h>
/*
  FFTW 3 behind the FFTW 2 names strack and strackw were written against.

  The trackers keep FFTW 2's struct-style complex type, with .re/.im accessed at a
  few hundred sites, so that type is kept here and cast to fftwf_complex at the
  handful of plan/execute calls in strackFFT.c.  The two are layout-identical: two
  consecutive floats.  fftw3.h is deliberately NOT included by this header, because
  it defines its own (double) fftw_complex which would clash; strackFFT.c is the
  only file that includes it.

  Plans are shared across threads and executed with fftwf_execute_dft(), which is
  thread safe on distinct, equally aligned arrays.  Plan creation is not, and is
  serialised inside strackPlan2d().  Every array that is passed to strackExec2d()
  must come from strackMallocComplex(), so its alignment matches what the plans
  were made on - a plain malloc'd array is not guaranteed to.

  Plans are looked up in FFTW wisdom, seeded once per host with FFTW_PATIENT and
  replayed on every later run; see strackFFT.c for the file location and why.
*/
typedef struct
{
	float re;
	float im;
} fftw_complex; /* same layout as fftwf_complex */
typedef float fftw_real;
typedef struct fftwf_plan_s *fftwnd_plan; /* an fftwf_plan, under the old name */

#ifndef FFTW_FORWARD
#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)
#endif

/* 2-D complex plan of nA rows by nR columns, out of place, in the given direction */
fftwnd_plan strackPlan2d(int32_t nA, int32_t nR, int32_t direction);
/* execute plan on in -> out; both must be from strackMallocComplex, contiguous, nA*nR long */
void strackExec2d(fftwnd_plan plan, void *in, void *out);
void strackDestroyPlan(fftwnd_plan plan);
/* n complex values, aligned as the plans require */
void *strackMallocComplex(size_t n);
void strackFreeComplex(void *p);
#endif
