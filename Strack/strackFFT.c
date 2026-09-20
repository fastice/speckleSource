/*
  strackFFT.c - FFTW 3 plans and execution for strack and strackw, with wisdom.

  This is the only file in the trackers that includes fftw3.h.  The trackers keep
  FFTW 2's struct-style complex type (see strackFFT.h); it is layout-identical to
  fftwf_complex and is cast here at the boundary.

  Planning rigor and wisdom, the same scheme as lstrack (landsatSource64):

  FFTW_ESTIMATE picks poor plans for the awkward sizes the trackers use (72, 96, 288
  are 9*8, 3*32, 9*32) - measured 1.7-1.8x slower than a measured plan at 72 and 96.
  But a measured plan is chosen by timing candidates at run time, so on a loaded
  machine two runs can pick differently and the offsets then differ in the last
  bit.  Wisdom settles both: search once, store the choice, replay it every later
  run - full speed and byte-reproducible output.

  The search runs once per host and whatever it finds is frozen in the file, so it
  searches as hard as it reasonably can (FFTW_PATIENT, a few seconds per size).  The
  file is $STRACK_WISDOM, else $HOME/.strackWisdom.<hostname>.  Per host because
  $HOME is shared between machines here: with one shared file, a host that could
  not replay the other's plan would re-measure and overwrite it, and the two would
  clobber each other on every run.  The file self-extends: a size not yet covered
  (new wr/wra in a par file) is measured and merged in, existing entries kept.
  Writing goes through a temp file and rename so concurrent tracker processes
  cannot tear each other's file.

  Seed it on a quiet machine with a single process before a bulk campaign; a
  measurement taken inside an 8-way batch is taken under contention, and the last
  writer wins.  Deleting the file breaks nothing - the next run re-seeds.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fftw3.h>

/* the same tokens strackFFT.h hands the trackers; fftw3.h agrees on the values */
#if FFTW_FORWARD != -1 || FFTW_BACKWARD != 1
#error "FFTW direction constants differ from the values strackFFT.h assumes"
#endif

static int32_t wisdomImported = 0; /* import once per process */
static int32_t wisdomDirty = 0;	   /* a plan was searched for; file needs writing */

static const char *strackWisdomPath(char *buf, size_t n)
{
	const char *envPath, *home;
	char hostName[256];
	envPath = getenv("STRACK_WISDOM");
	if (envPath != NULL && strlen(envPath) > 0)
	{
		return (envPath);
	}
	home = getenv("HOME");
	if (home == NULL)
	{
		return (NULL);
	}
	if (gethostname(hostName, sizeof(hostName)) != 0)
	{
		strcpy(hostName, "unknown");
	}
	hostName[sizeof(hostName) - 1] = '\0';
	snprintf(buf, n, "%s/.strackWisdom.%s", home, hostName);
	return (buf);
}

static void saveStrackWisdom(const char *path)
{
	char tmpPath[2176]; /* room for the path plus ".tmp.<pid>" */
	if (path == NULL)
	{
		return;
	}
	snprintf(tmpPath, sizeof(tmpPath), "%s.tmp.%i", path, (int)getpid());
	if (fftwf_export_wisdom_to_filename(tmpPath) == 0)
	{
		fprintf(stderr, "warning: could not write FFTW wisdom to %s\n", tmpPath);
		return;
	}
	if (rename(tmpPath, path) != 0)
	{
		fprintf(stderr, "warning: could not install FFTW wisdom at %s\n", path);
		remove(tmpPath);
		return;
	}
	fprintf(stderr, "wrote FFTW wisdom to %s\n", path);
}

/*
  2-D complex out-of-place plan, nA rows by nR columns.  Serialised: the FFTW
  planner holds global state and is not thread safe, and the trackers create their
  plans from inside an omp parallel region.
*/
struct fftwf_plan_s *strackPlan2d(int32_t nA, int32_t nR, int32_t direction)
{
	struct fftwf_plan_s *plan = NULL;
	fftwf_complex *in, *out;
	const char *wisdomPath;
	char wisdomBuf[2048];
	size_t n = (size_t)nA * (size_t)nR;
#pragma omp critical(strackFftPlanner)
	{
		wisdomPath = strackWisdomPath(wisdomBuf, sizeof(wisdomBuf));
		if (wisdomImported == 0)
		{
			if (wisdomPath != NULL)
			{
				fftwf_import_wisdom_from_filename(wisdomPath);
			}
			wisdomImported = 1;
		}
		/* scratch arrays for planning: the plans are only ever executed through
		   fftwf_execute_dft on the trackers' own strackMallocComplex arrays, which
		   share this alignment, so these can go once the plan exists */
		in = fftwf_malloc(sizeof(fftwf_complex) * n);
		out = fftwf_malloc(sizeof(fftwf_complex) * n);
		if (in == NULL || out == NULL)
		{
			fprintf(stderr, "strackPlan2d: could not allocate planning arrays for %i x %i\n", nA, nR);
			exit(1);
		}
		/* replay the stored plan if the wisdom covers this transform */
		plan = fftwf_plan_dft_2d(nA, nR, in, out, direction, FFTW_PATIENT | FFTW_WISDOM_ONLY);
		if (plan == NULL)
		{
			fprintf(stderr, "\033[33mno FFTW wisdom for %i x %i at %s - searching now (once per host); "
							"this run may differ in the last bit from later runs that replay it\033[0m\n",
					nA, nR, wisdomPath == NULL ? "(no path)" : wisdomPath);
			plan = fftwf_plan_dft_2d(nA, nR, in, out, direction, FFTW_PATIENT);
			if (plan == NULL)
			{
				fprintf(stderr, "strackPlan2d: could not create a %i x %i plan\n", nA, nR);
				exit(1);
			}
			wisdomDirty = 1;
			saveStrackWisdom(wisdomPath); /* merged export: earlier sizes are kept */
		}
		fftwf_free(in);
		fftwf_free(out);
	}
	return (plan);
}

void strackExec2d(struct fftwf_plan_s *plan, void *in, void *out)
{
	fftwf_execute_dft(plan, (fftwf_complex *)in, (fftwf_complex *)out);
}

void strackDestroyPlan(struct fftwf_plan_s *plan)
{
	if (plan == NULL)
	{
		return;
	}
#pragma omp critical(strackFftPlanner)
	{
		fftwf_destroy_plan(plan);
	}
}

void *strackMallocComplex(size_t n)
{
	return (fftwf_malloc(sizeof(fftwf_complex) * n));
}

void strackFreeComplex(void *p)
{
	if (p != NULL)
	{
		fftwf_free(p);
	}
}
