CC =		gcc
CFLAGS  =	$(FLAGS) -c  -I$(INCLUDEPATH) -D_LARGEFILE64_SOURCE -D_FILE64 -D_LARGEFILE_SOURCE -DHPLARGEFILES

OBJS= strack.o parseTrack.o parseBase.o parseInitialOffsets.o speckleTrack.o \
      sTrackOut.o getInt.o getMask.o parsePar.o writeVrt.o readBothOffsetsStrack.o

MACHTYPE = $(shell uname -m)
OSTYPE = $(shelf uname -s)

$(shell mkdir -p $(MACHTYPE)-$(OSTYPE))


.KEEP_STATE:
all:	$(OBJS)
	mv *.o $(MACHTYPE)-$(OSTYPE)

