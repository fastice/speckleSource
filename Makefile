C =		gcc
ROOTDIR =	/Users/ian
PROGDIR =       $(ROOTDIR)/progs/GIT
INCLUDEPATH =	$(ROOTDIR)/progs/GIT
BINDIR =	$(IHOME)/bin/$(MACHTYPE)
#
CFLAGS =	'-O3 -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS) '
CCFLAGS =  '-O3 -m32 -D$(MACHTYPE) $(COMPILEFLAGS) '
#-Wunused-variable'

CCFLAGS1= '-O3 -march=native'
# uncomment to debug
#CFLAGS =	'-g -m32 -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g -m32 -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'



FFLAGS  =	 -g
ERS1CODEDIR =	ers1Code

COMMON=	$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/earthRadiusFunctions.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/parseInputFile.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/lltoxy1.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readOldPar.o




FFT =	$(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_1.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_2.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_3.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_4.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_5.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_6.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_7.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_8.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_9.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_10.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_11.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_12.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_13.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_14.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_15.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fn_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_2.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_3.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_4.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_5.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_6.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_7.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_8.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_9.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_10.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftw_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_1.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_2.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_3.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_4.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_5.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_6.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_7.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_8.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_9.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_10.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_11.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_12.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_13.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_14.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_15.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fni_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_2.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_3.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_4.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_5.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_6.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_7.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_8.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_9.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_10.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/ftwi_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/timer.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/config.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/planner.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/twiddle.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/fftwnd.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/wisdom.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/wisdomio.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/putils.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/rader.o  $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/malloc.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/generic.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(MACHTYPE)-$(OSTYPE)/executor.o

RDF =		$(PROGDIR)/rdfSource/rdfRoutines/$(MACHTYPE)-$(OSTYPE)/SRTMrdf.o

STANDARD =	$(PROGDIR)/clib/$(MACHTYPE)-$(OSTYPE)/standard.o

RECIPES  =	$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/polint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/nrutil.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/ratint.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/four1.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdfit.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdcmp.o \
		$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svdvar.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/svbksb.o $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/pythag.o \
                $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/hunt.o  $(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/gaussj.o \
		$(PROGDIR)/cRecipes/$(MACHTYPE)-$(OSTYPE)/covsrt.o


UNWRAP = $(PROGDIR)/unwrapSource/unWrap/$(MACHTYPE)-$(OSTYPE)/labelRegions.o



TARGETS = strack  strackw  cullst cullls

all: $(TARGETS)

#******************************************************************************************************************
#*********************************************strack **************************************************************
#******************************************************************************************************************

STRACK  =	Strack/$(MACHTYPE)-$(OSTYPE)/parseTrack.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseInitialOffsets.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/speckleTrack.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/sTrackOut.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/getInt.o \
	 	Strack/$(MACHTYPE)-$(OSTYPE)/getMask.o \
                Strack/$(MACHTYPE)-$(OSTYPE)/parseBase.o \
                Strack/$(MACHTYPE)-$(OSTYPE)/parsePar.o

STRACKDIRS =	Strack 

strack:
	@for i in ${STRACKDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc -m32   $(CCFLAGS1) \
                Strack/$(MACHTYPE)-$(OSTYPE)/strack.o $(STRACK)  $(STANDARD) $(RECIPES)  $(RDF) $(FFT) $(COMMON)  \
                -lm  -o $(BINDIR)/strack

#******************************************************************************************************************
#*********************************************strackw **************************************************************
#******************************************************************************************************************

STRACKW	=	Strackw/$(MACHTYPE)-$(OSTYPE)/corrTrackFast.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseTrack.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseInitialOffsets.o \
	 	Strack/$(MACHTYPE)-$(OSTYPE)/getMask.o

STRACKWDIRS =	Strackw Strack

strackw:
	@for i in ${STRACKWDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc -m32   $(CCFLAGS1) \
		Strackw/$(MACHTYPE)-$(OSTYPE)/strackw.o $(STRACKW) $(STANDARD) $(RECIPES) $(RDF) $(FFT)  $(COMMON) \
		Strack/$(MACHTYPE)-$(OSTYPE)/parsePar.o \
                -lm  -o $(BINDIR)/strackw


#******************************************************************************************************************
#********************************************cullst **************************************************************
#******************************************************************************************************************

CULLST  =	Cullst/$(MACHTYPE)-$(OSTYPE)/loadCullData.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullSTData.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullStats.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullSmooth.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullIslands.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/writeCullData.o

CULLSTDIRS =	Cullst

cullst:
	@for i in ${CULLSTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc -m32 $(CCFLAGS1) \
                Cullst/$(MACHTYPE)-$(OSTYPE)/cullst.o $(CULLST) $(STANDARD) $(RECIPES) $(UNWRAP) \
                -lm  -o $(BINDIR)/cullst


#******************************************************************************************************************
#********************************************cullst **************************************************************
#******************************************************************************************************************

CULLLS  =	Cullls/$(MACHTYPE)-$(OSTYPE)/loadLSCullData.o \
			Cullls/$(MACHTYPE)-$(OSTYPE)/cullLSData.o \
			Cullls/$(MACHTYPE)-$(OSTYPE)/cullLSStats.o \
			Cullls/$(MACHTYPE)-$(OSTYPE)/cullLSSmooth.o \
			Cullls/$(MACHTYPE)-$(OSTYPE)/cullLSIslands.o

LANDSATMOSAIC =		$(PROGDIR)/mosaicSource/landsatMosaic/$(MACHTYPE)-$(OSTYPE)/readLSOffsets.o
CULLLSDIRS =	Cullls

cullls:
	@for i in ${CULLLSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/speckleSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                Cullls/$(MACHTYPE)-$(OSTYPE)/cullls.o $(CULLLS) $(STANDARD) $(RECIPES) $(UNWRAP) $(LANDSATMOSAIC) \
                -lm  -o $(BINDIR)/cullls

TESTDIRS =	test
testg:
	@for i in ${TESTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/speckleSource; \
		); done
		gcc -m32 $(CCFLAGS1) \
                test/$(MACHTYPE)-$(OSTYPE)/test.o $(STANDARD) $(RECIPES)  \
                -lm  -o $(BINDIR)/testg


