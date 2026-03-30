
#include "strack.h"
#include "math.h"
#include <sys/types.h>
static double computeSingleDop(fftw_complex **image, int32_t nAz, int32_t nRg);
/*
  Estimate Doppler centroid (Hz) from RSLC block.

  image1: complex image addressed as image1[j][i]
          j = 0..nr-1 (range/cols)
          i = 0..na-1 (az/rows)

  i0, j0: starting indices (row, col)
  nAz:    number of azimuth rows to use (e.g., 64)
  nRg:    number of range cols to average (e.g., nRangeBins)
  prf_hz: PRF (Hz)

  Returns: Doppler centroid in Hz. If window invalid, returns NAN.
*/
static int iround(double x)
{
    return (x >= 0.0) ? (int)(x + 0.5) : (int)(x - 0.5);
}

void estDopCarrier1(TrackParams *trackPar,  fftw_complex **image1, fftw_complex **image2, int32_t nAz, int32_t nRg)
{
    int i, j;
	double angle1, angle2;
	int32_t shift;
  	angle1 = computeSingleDop(image1, nAz, nRg);
	angle2 = computeSingleDop(image2, nAz, nRg);
	//fprintf(stderr, "\nxxx %f %f %f\n", angle1, angle2, angle1-angle2);
	shift = iround(nAz/2 *(angle1 + angle2)*.5/PI);
	if(shift < 0) shift += nAz;
	trackPar->azShift = shift;
	// Normalize back to complex window size to be consistent with cmp case
	// Will get renormalized for amplitude case
	trackPar->azShift = (int)((float)trackPar->azShift * (float)trackPar->wA/(float)nAz);
	//printf(stdout, "\n+ %lf %lf %i %i\n", angle1, angle1/PI, nAz, trackPar->azShift);
}

static double computeSingleDop(fftw_complex **image, int32_t nAz, int32_t nRg)
{
	int32_t i, j;
	int32_t i0 = 0, j0 = 0;
	double a_re, a_im, b_re, b_im;
   	double mean_re=0., mean_im=0.;
	double sum_re = 0.0;
	double sum_im = 0.0;
 	long double denom = (long double)nRg * (long double)nAz;
	double z_re, z_im;
	double angle_rad;
	// Optional but often helpful: remove a single complex mean over the whole window
	// to reduce DC/leakage. (If you don't want this, delete this block and set mean=0.)
	for (j = j0; j < j0 + nRg; ++j) {
        for (i = i0; i < i0 + nAz; ++i) {
            fftw_complex *p1 = &(image[i][j]); /* image1[j][i] flattened */
			fftw_complex *p2 = &(image[i][j]); /* image1[j][i] flattened */
            mean_re += (long double)p1->re;
            mean_im += (long double)p1->im;
        }
    }
    mean_re /= denom;
    mean_im /= denom;
	// Now dop
	for (j = j0; j < j0 + nRg; ++j) {
		for (i = i0; i < i0 + nAz - 1; ++i) {
			fftw_complex *p0 = &(image[i][j]);
			fftw_complex *p1 = &(image[i + 1][j]);

			/* s0 = p0 - mean, s1 = p1 - mean */
			a_re = (long double)p1->re - mean_re;
			a_im = (long double)p1->im - mean_im;
			b_re = (long double)p0->re - mean_re;
			b_im = (long double)p0->im - mean_im;

			/* a * conj(b): re = ar*br + ai*bi ; im = ai*br - ar*bi */
			z_re = a_re * b_re + a_im * b_im;
			z_im = a_im * b_re - a_re * b_im;

			sum_re += z_re;
			sum_im += z_im;
		}
	}
	if (sum_re == 0.0L && sum_im == 0.0L) return (0.0/0.0); /* NAN */
	angle_rad = atan2((double)sum_im, (double)sum_re);
	return angle_rad;
}
