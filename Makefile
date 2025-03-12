#
# Modify this section to point to where stuff is.
# Current names are for a specific file system.
# ROOTDIR: root directory for code (e.g. /home/username directory).
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
$(info  "USER is $(USER)")
MACHTYPE = $(shell uname -m)
OSTYPE = $(shell uname -s)
$(info "MACHTYPE/OSTYPE $(MACHTYPE) $(OSTYPE)")
#
# Default rootdir
ifneq ($(ROOTDIR)),)
	ROOTDIR =	$(dir $(CURDIR))
endif
$(info ROOTDIR="$(ROOTDIR)")
#
# Default root for source code
ifneq ($(PROGDIR)),)
	PROGDIR =       $(dir $(CURDIR))
endif
$(info PROGDIR ="$(PROGDIR)")
#
# Default location root for compiled programs
ifneq ($(BINHOME)),)
	BINHOME =		~$(USER)
endif
# For historical reasons, can compile with 32-bit memory model using MEM=-m32
# In almost all cases, should be compiled as 64bit.
ifneq ($(MEM), -m32)
	BINNAME=	$(MACHTYPE)
	FFTDIR = $(MACHTYPE)-$(OSTYPE)
else
	BINNAME =	i386
	FFTDIR = i386-$(OSTYPE)
endif
$(info "FFTDIR = $(FFTDIR)")
FFTHOME = 	$(PROGDIR)/fft/fftw-2.1.5/fftw
#
# Default binary directory
ifneq ($(BINDIR)),)
	BINDIR =	$(BINHOME)/bin/$(BINNAME)
endif
$(info "BINDIR = $(BINDIR)")
#
# Create bin dir if it doesn't exist
$(shell mkdir -p $(BINDIR))
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
#
CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
CCFLAGS =  '-O3 $(MEM) $(COMPILEFLAGS) '
GDAL = -lgdal -lcurl  -lsqlite3 -llzma -lpoppler -lopenjp2 -lssh2 -llcms2
#
CCFLAGS1= -O3 
#-no-pie
# uncomment to debug
#CFLAGS =	'-g $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
#CCFLAGS =  '-g $(MEM) -D$(MACHTYPE) $(COMPILEFLAGS)'
#CCFLAGS1= '-g'
#
ifneq ("$(OSTYPE)", "Darwin")
	NOPIE =	-no-pie
	GDAL = -lgdal -lcurl  -lsqlite3 -llzma -lpoppler -lopenjp2 -lssh2 -llcms2
	CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS)'
	CCFLAGS =  '-O3 $(MEM) $(COMPILEFLAGS) '
else
	GDALLIB = /opt/homebrew/lib
	GDALINCLUDE = /opt/homebrew/include
	GDAL = -lgdal -L/opt/homebrew/lib
	CFLAGS =	'-O3 $(MEM) -I$(INCLUDEPATH) $(COMPILEFLAGS) -I$(GDALINCLUDE)'
	CCFLAGS =  '-O3 $(MEM) $(COMPILEFLAGS) -I$(GDALINCLUDE)'
endif
#
# ******** SHOULD NOT NEED TO MODIFY BELOW HERE *********
#

FFLAGS  =	 -g
ERS1CODEDIR =	ers1Code

COMMON=	$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/addIrregData.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/bilinearInterp.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/computeHeading.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/computePhiZ.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/computeScale.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/computeTiePoints.o \
	    	$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/computeXYangle.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/earthRadiusFunctions.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/geojsonCode.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getDataStringSpecial.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getBaseline.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getHeight.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getMVhInputFile.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getRegion.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getShelfMask.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/getXYHeight.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/groundRangeToLLNew.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/initMatrix.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/initRoutines.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/interpOffsets.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/interpPhaseImage.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/interpTideDiff.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/interpVCorrect.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/interpXYDEM.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/julianDay.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/llToImageNew.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/lltoxy.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/lltoxy1.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/outputGeocodedImage.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/parseInputFile.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/parseIrregFile.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/polintVec.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/rangeAzimuthToLL.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readOffsets.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readOldPar.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readShelf.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readTiePoints.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/readXYDEM.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/rotateFlowDirectionToRA.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/rotateFlowDirectionToXY.o \
	        $(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/scalingFunctions.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/smlocateZD.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/svBase.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/vectorFunc.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/xyGetZandSlope.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/xytoll1.o \
			$(PROGDIR)/mosaicSource/common/$(MACHTYPE)-$(OSTYPE)/xytoll.o

FFT = $(FFTHOME)/$(FFTDIR)/fn_1.o $(FFTHOME)/$(FFTDIR)/fn_2.o \
       $(FFTHOME)/$(FFTDIR)/fn_3.o $(FFTHOME)/$(FFTDIR)/fn_4.o \
       $(FFTHOME)/$(FFTDIR)/fn_5.o $(FFTHOME)/$(FFTDIR)/fn_6.o \
       $(FFTHOME)/$(FFTDIR)/fn_7.o $(FFTHOME)/$(FFTDIR)/fn_8.o \
       $(FFTHOME)/$(FFTDIR)/fn_9.o $(FFTHOME)/$(FFTDIR)/fn_10.o \
       $(FFTHOME)/$(FFTDIR)/fn_11.o $(FFTHOME)/$(FFTDIR)/fn_12.o \
       $(FFTHOME)/$(FFTDIR)/fn_13.o $(FFTHOME)/$(FFTDIR)/fn_14.o \
       $(FFTHOME)/$(FFTDIR)/fn_15.o $(FFTHOME)/$(FFTDIR)/fn_16.o \
       $(FFTHOME)/$(FFTDIR)/fn_32.o $(FFTHOME)/$(FFTDIR)/fn_64.o \
       $(FFTHOME)/$(FFTDIR)/ftw_2.o $(FFTHOME)/$(FFTDIR)/ftw_3.o \
       $(FFTHOME)/$(FFTDIR)/ftw_4.o $(FFTHOME)/$(FFTDIR)/ftw_5.o \
       $(FFTHOME)/$(FFTDIR)/ftw_6.o $(FFTHOME)/$(FFTDIR)/ftw_7.o \
       $(FFTHOME)/$(FFTDIR)/ftw_8.o $(FFTHOME)/$(FFTDIR)/ftw_9.o \
       $(FFTHOME)/$(FFTDIR)/ftw_10.o $(FFTHOME)/$(FFTDIR)/ftw_16.o \
       $(FFTHOME)/$(FFTDIR)/ftw_32.o $(FFTHOME)/$(FFTDIR)/ftw_64.o \
       $(FFTHOME)/$(FFTDIR)/fni_1.o $(FFTHOME)/$(FFTDIR)/fni_2.o \
       $(FFTHOME)/$(FFTDIR)/fni_3.o $(FFTHOME)/$(FFTDIR)/fni_4.o \
       $(FFTHOME)/$(FFTDIR)/fni_5.o $(FFTHOME)/$(FFTDIR)/fni_6.o \
       $(FFTHOME)/$(FFTDIR)/fni_7.o $(FFTHOME)/$(FFTDIR)/fni_8.o \
       $(FFTHOME)/$(FFTDIR)/fni_9.o $(FFTHOME)/$(FFTDIR)/fni_10.o \
       $(FFTHOME)/$(FFTDIR)/fni_11.o $(FFTHOME)/$(FFTDIR)/fni_12.o \
       $(FFTHOME)/$(FFTDIR)/fni_13.o $(FFTHOME)/$(FFTDIR)/fni_14.o \
       $(FFTHOME)/$(FFTDIR)/fni_15.o $(FFTHOME)/$(FFTDIR)/fni_16.o \
       $(FFTHOME)/$(FFTDIR)/fni_32.o $(FFTHOME)/$(FFTDIR)/fni_64.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_2.o $(FFTHOME)/$(FFTDIR)/ftwi_3.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_4.o $(FFTHOME)/$(FFTDIR)/ftwi_5.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_6.o $(FFTHOME)/$(FFTDIR)/ftwi_7.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_8.o $(FFTHOME)/$(FFTDIR)/ftwi_9.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_10.o $(FFTHOME)/$(FFTDIR)/ftwi_16.o \
       $(FFTHOME)/$(FFTDIR)/ftwi_32.o $(FFTHOME)/$(FFTDIR)/ftwi_64.o \
       $(FFTHOME)/$(FFTDIR)/timer.o $(FFTHOME)/$(FFTDIR)/config.o \
       $(FFTHOME)/$(FFTDIR)/planner.o $(FFTHOME)/$(FFTDIR)/twiddle.o \
       $(FFTHOME)/$(FFTDIR)/fftwnd.o $(FFTHOME)/$(FFTDIR)/wisdom.o \
       $(FFTHOME)/$(FFTDIR)/wisdomio.o $(FFTHOME)/$(FFTDIR)/putils.o \
       $(FFTHOME)/$(FFTDIR)/rader.o $(FFTHOME)/$(FFTDIR)/malloc.o \
       $(FFTHOME)/$(FFTDIR)/generic.o $(FFTHOME)/$(FFTDIR)/executor.o

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

GDALIO = 	$(PROGDIR)/gdalIO/gdalIO/$(MACHTYPE)-$(OSTYPE)/gdalIO.o \
			$(PROGDIR)/gdalIO/gdalIO/$(MACHTYPE)-$(OSTYPE)/dictionaryCode.o \
			$(PROGDIR)/gdalIO/gdalIO/$(MACHTYPE)-$(OSTYPE)/tiffWriteCode.o


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
		Strack/$(MACHTYPE)-$(OSTYPE)/writeVrt.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/readBothOffsetsStrack.o \
        Strack/$(MACHTYPE)-$(OSTYPE)/parsePar.o

STRACKDIRS =	Strack $(PROGDIR)/rdfSource/rdfRoutines $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/mosaicSource/common $(FFTHOME) \
	$(PROGDIR)/rdfSource/rdfRoutines

strack:
	@for i in ${STRACKDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		g++ $(MEM) $(CCFLAGS1) $(NOPIE) \
                Strack/$(MACHTYPE)-$(OSTYPE)/strack.o $(STRACK)  $(STANDARD) $(RECIPES)  $(RDF) $(FFT) $(COMMON) $(GDALIO) \
                -lm $(LDFLAGS)  $(GDAL) -o $(BINDIR)/strack

#******************************************************************************************************************
#*********************************************strackw **************************************************************
#******************************************************************************************************************

STRACKW	=	Strackw/$(MACHTYPE)-$(OSTYPE)/corrTrackFast.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseTrack.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/parseInitialOffsets.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/sTrackOut.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/writeVrt.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/readBothOffsetsStrack.o \
	 	Strack/$(MACHTYPE)-$(OSTYPE)/getMask.o

STRACKWDIRS =	Strackw Strack $(PROGDIR)/rdfSource/rdfRoutines $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/clib $(PROGDIR)/cRecipes $(PROGDIR)/mosaicSource/common $(FFTHOME) \
	$(PROGDIR)/rdfSource/rdfRoutines

strackw:
	@for i in ${STRACKWDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		g++ $(MEM)   $(CCFLAGS1) $(NOPIE) \
		Strackw/$(MACHTYPE)-$(OSTYPE)/strackw.o $(STRACKW) $(STANDARD) $(RECIPES) $(RDF) $(FFT)  $(COMMON) $(GDALIO) \
		Strack/$(MACHTYPE)-$(OSTYPE)/parsePar.o \
                -lm $(LDFLAGS)  $(GDAL) -o $(BINDIR)/strackw


#******************************************************************************************************************
#********************************************cullst **************************************************************
#******************************************************************************************************************

CULLST  =	Cullst/$(MACHTYPE)-$(OSTYPE)/loadCullData.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullSTData.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullStats.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullSmooth.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/cullIslands.o \
		Strack/$(MACHTYPE)-$(OSTYPE)/writeVrt.o \
		Cullst/$(MACHTYPE)-$(OSTYPE)/writeCullData.o

CULLSTDIRS =	Cullst $(PROGDIR)/clib $(PROGDIR)/gdalIO/gdalIO $(PROGDIR)/cRecipes  $(PROGDIR)/unwrapSource/unWrap

cullst:
	@for i in ${CULLSTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR); \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                Cullst/$(MACHTYPE)-$(OSTYPE)/cullst.o $(CULLST) $(STANDARD) $(RECIPES) $(COMMON) $(UNWRAP) $(GDALIO)\
                -lm  $(LDFLAGS) $(GDAL)  -o $(BINDIR)/cullst


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
		g++ $(MEM) $(CCFLAGS1) \
                Cullls/$(MACHTYPE)-$(OSTYPE)/cullls.o $(CULLLS) $(STANDARD) $(RECIPES) $(UNWRAP) $(LANDSATMOSAIC) \
                -lm  $(GDAL) -o $(BINDIR)/cullls

TESTDIRS =	test
testg:
	@for i in ${TESTDIRS}; do \
		( 	echo "<<< Descending in directory: $$i >>>"; \
	                cd $$i; \
			make FLAGS=$(CCFLAGS) INCLUDEPATH=$(INCLUDEPATH) PAF=0;  \
			cd $(PROGDIR)/speckleSource; \
		); done
		g++ $(MEM) $(CCFLAGS1) \
                test/$(MACHTYPE)-$(OSTYPE)/test.o $(STANDARD) $(RECIPES)  \
                -lm  -o $(BINDIR)/testg


