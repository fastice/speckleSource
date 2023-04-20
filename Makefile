#
# Modify this section to point to where stuff is.
# Current names are for a specific file system.
# ROOTDIR: root directory for code (e.g. /Users/username directory). Likely should change for linux
# PROGDIR: location for top source code directory (default ROOTDIR/progs/GIT64)
# BINHOME: root directory for binaries
# BINNAME: archetecture dependent basename for bin dir
# BINDIR: directory for binaries (default BINHOME/bin/BINNAME) (will create if doesn't exist)
# INCLUDEPATH: include path (default PROGDIR anything else could cause a problem)
# Various directors can be overridden with environment variable or from make command
# make BINHOME=/a/path/to/binhome
#
# Base directory for code
USER =	$(shell id -u -n)
#
# Default rootdir
ifneq ($(ROOTDIR)),)
	ROOTDIR =	/Users/$(USER)
endif
$(info ROOTDIR="$(ROOTDIR)")
# Default root for source code
ifneq ($(PROGDIR)),)
	PROGDIR =       $(ROOTDIR)/progs/GIT64
endif
$(info PROGDIR ="$(PROGDIR)")
#
# Default location root for compiled programs
ifneq ($(BINHOME)),)
	BINHOME =		~$(USER)
endif
$(info BINHOME="$(BINHOME)")
#
# For historical reasons, can compile with 32-bit memory model using MEM=-m32
# In almost all cases, should be compiled as 64bit.
ifneq ($(MEM),-m32)
	BINNAME=	$(MACHTYPE)
	FFTDIR = $(MACHTYPE)-$(OSTYPE)
else
	BINNAME =	i386
	FFTDIR = i386-$(OSTYPE)
endif
$(info BINNAME="$(BINNAME)")
#
# Default binary directory
ifneq ($(BINDIR)),)
	BINDIR =	$(BINHOME)/bin/$(BINNAME)
endif
$(info BINDIR="$(BINDIR)")
#
# Create bin dir if it doesn't exist
$(shell mkdir -p $(BINDIR))
#
#
# Default include path
ifneq ($(INCLUDEPATH)),)
	INCLUDEPATH =	$(PROGDIR)
endif
$(info INCLUDEPATH ="$(INCLUDEPATH)")
#
# Compiler stuff
#
C =		gcc
CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 $(MEM) -D$(MACHTYPE) $(COMPILEFLAGS)'
#-Wunused-variable'

CCFLAGS1= '-O3'
# uncomment to debug
#CFLAGS =	'-g $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g $(MEM) -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'
#
# ******** SHOULD NOT NEED TO MODIFY BELOW HERE *********
#

FFLAGS  =	 -g
ERS1CODEDIR =	ers1Code

COMMON=	$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/earthRadiusFunctions.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/parseInputFile.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/lltoxy1.o \
	                $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readOldPar.o


FFT =	$(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_1.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_2.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_3.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_4.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_5.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_6.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_7.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_8.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_9.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_10.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_11.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_12.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_13.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_14.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_15.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fn_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_2.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_3.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_4.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_5.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_6.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_7.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_8.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_9.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_10.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftw_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_1.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_2.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_3.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_4.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_5.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_6.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_7.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_8.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_9.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_10.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_11.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_12.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_13.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_14.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_15.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fni_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_2.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_3.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_4.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_5.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_6.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_7.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_8.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_9.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_10.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_16.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_32.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/ftwi_64.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/timer.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/config.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/planner.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/twiddle.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/fftwnd.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/wisdom.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/wisdomio.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/putils.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/rader.o  $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/malloc.o \
        $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/generic.o $(PROGDIR)/fft/fftw-2.1.5/fftw/$(FFTDIR)/executor.o

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

STRACKDIRS =	Strack $(PROGDIR)/rdfSource/rdfRoutines $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/mosaicSource/common

strack:
	@for i in ${STRACKDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc $(MEM) $(CCFLAGS1) -no-pie \
                Strack/$(MACHTYPE)-$(OSTYPE)/strack.o $(STRACK)  $(STANDARD) $(RECIPES)  $(RDF) $(FFT) $(COMMON)  \
                -lm  -o $(BINDIR)/strack

#******************************************************************************************************************
#*********************************************strackw **************************************************************
#******************************************************************************************************************

STRACKW	=	Strackw/$(MACHTYPE)-$(OSTYPE)/corrTrackFast.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseTrack.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseInitialOffsets.o \
	 	Strack/$(MACHTYPE)-$(OSTYPE)/getMask.o

STRACKWDIRS =	Strackw Strack $(PROGDIR)/rdfSource/rdfRoutines $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/mosaicSource/common

strackw:
	@for i in ${STRACKWDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc $(MEM)   $(CCFLAGS1) -no-pie \
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

CULLSTDIRS =	Cullst $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/unwrapSource/unWrap

cullst:
	@for i in ${CULLSTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		gcc $(MEM) $(CCFLAGS1) \
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
CULLLSDIRS =	Cullls $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/unwrapSource/unWrap $(PROGDIR)/mosaicSource/landsatMosaic

cullls:
	@for i in ${CULLLSDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/speckleSource; \
		); done
		gcc $(MEM) $(CCFLAGS1) \
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
		gcc $(MEM) $(CCFLAGS1) \
                test/$(MACHTYPE)-$(OSTYPE)/test.o $(STANDARD) $(RECIPES)  \
                -lm  -o $(BINDIR)/testg


